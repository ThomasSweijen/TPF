// 2009 © Václav Šmilauer <eudoxos@arcig.cz>


#include<lib/base/Math.hpp>
#include<lib/base/openmp-accu.hpp>

#include<core/Engine.hpp>

#include<pkg/common/Dispatching.hpp>
#include<pkg/common/Callbacks.hpp>
#include<pkg/dem/SpherePack.hpp>
#include<pkg/common/KinematicEngines.hpp>
#ifdef YADE_OPENGL
	#include<pkg/common/GLDrawFunctors.hpp>
	#include<pkg/common/OpenGLRenderer.hpp>
#endif
#include<pkg/common/MatchMaker.hpp>

// move this to the miniEigen wrapper later

/* two-way se3 handling */
struct custom_se3_to_tuple{
	static PyObject* convert(const Se3r& se3){
		boost::python::tuple ret=boost::python::make_tuple(se3.position,se3.orientation);
		return boost::python::incref(ret.ptr());
	}
};
struct custom_Se3r_from_seq{
	custom_Se3r_from_seq(){
		 boost::python::converter::registry::push_back(&convertible,&construct,boost::python::type_id<Se3r>());
	}
	static void* convertible(PyObject* obj_ptr){
		 if(!PySequence_Check(obj_ptr)) return 0;
		 if(PySequence_Size(obj_ptr)!=2 && PySequence_Size(obj_ptr)!=7) return 0;
		 return obj_ptr;
	}
	static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data){
		void* storage=((boost::python::converter::rvalue_from_python_storage<Se3r>*)(data))->storage.bytes;
		new (storage) Se3r; Se3r* se3=(Se3r*)storage;
		if(PySequence_Size(obj_ptr)==2){ // from vector and quaternion
			se3->position=boost::python::extract<Vector3r>(PySequence_GetItem(obj_ptr,0));
			se3->orientation=boost::python::extract<Quaternionr>(PySequence_GetItem(obj_ptr,1));
		} else if(PySequence_Size(obj_ptr)==7){ // 3 vector components, 3 axis components, angle
			se3->position=Vector3r(boost::python::extract<Real>(PySequence_GetItem(obj_ptr,0)),boost::python::extract<Real>(PySequence_GetItem(obj_ptr,1)),boost::python::extract<Real>(PySequence_GetItem(obj_ptr,2)));
			Vector3r axis=Vector3r(boost::python::extract<Real>(PySequence_GetItem(obj_ptr,3)),boost::python::extract<Real>(PySequence_GetItem(obj_ptr,4)),boost::python::extract<Real>(PySequence_GetItem(obj_ptr,5)));
			Real angle=boost::python::extract<Real>(PySequence_GetItem(obj_ptr,6));
			se3->orientation=Quaternionr(AngleAxisr(angle,axis));
		} else throw std::logic_error(__FILE__ ": First, the sequence size for Se3r object was 2 or 7, but now is not? (programming error, please report!");
		data->convertible=storage;
	}
};


struct custom_OpenMPAccumulator_to_float{ static PyObject* convert(const OpenMPAccumulator<Real>& acc){ return boost::python::incref(PyFloat_FromDouble(acc.get())); } };
struct custom_OpenMPAccumulator_from_float{
	custom_OpenMPAccumulator_from_float(){  boost::python::converter::registry::push_back(&convertible,&construct,boost::python::type_id<OpenMPAccumulator<Real> >()); }
	static void* convertible(PyObject* obj_ptr){ return PyFloat_Check(obj_ptr) ? obj_ptr : 0; }
	static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data){ void* storage=((boost::python::converter::rvalue_from_python_storage<OpenMPAccumulator<Real> >*)(data))->storage.bytes; new (storage) OpenMPAccumulator<Real>; ((OpenMPAccumulator<Real>*)storage)->set(boost::python::extract<Real>(obj_ptr)); data->convertible=storage; }
};
struct custom_OpenMPAccumulator_to_int  { static PyObject* convert(const OpenMPAccumulator<int>& acc){ return boost::python::incref(PyInt_FromLong((long)acc.get())); } };
struct custom_OpenMPAccumulator_from_int{
	custom_OpenMPAccumulator_from_int(){  boost::python::converter::registry::push_back(&convertible,&construct,boost::python::type_id<OpenMPAccumulator<int> >()); }
	static void* convertible(PyObject* obj_ptr){ return PyInt_Check(obj_ptr) ? obj_ptr : 0; }
	static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data){ void* storage=((boost::python::converter::rvalue_from_python_storage<OpenMPAccumulator<int> >*)(data))->storage.bytes; new (storage) OpenMPAccumulator<int>; ((OpenMPAccumulator<int>*)storage)->set(boost::python::extract<int>(obj_ptr)); data->convertible=storage; }
};

template<typename T>
struct custom_vvector_to_list{
	static PyObject* convert(const std::vector<std::vector<T> >& vv){
		boost::python::list ret; FOREACH(const std::vector<T>& v, vv){
			boost::python::list ret2;
			FOREACH(const T& e, v) ret2.append(e);
			ret.append(ret2);
		}
		return boost::python::incref(ret.ptr());
	}
};

template<typename containedType>
struct custom_list_to_list{
	static PyObject* convert(const std::list<containedType>& v){
		boost::python::list ret; FOREACH(const containedType& e, v) ret.append(e);
		return boost::python::incref(ret.ptr());
	}
};
/*** c++-list to python-list */
template<typename containedType>
struct custom_vector_to_list{
	static PyObject* convert(const std::vector<containedType>& v){
		boost::python::list ret; FOREACH(const containedType& e, v) ret.append(e);
		return boost::python::incref(ret.ptr());
	}
};
template<typename containedType>
struct custom_vector_from_seq{
	custom_vector_from_seq(){ boost::python::converter::registry::push_back(&convertible,&construct,boost::python::type_id<std::vector<containedType> >()); }
	static void* convertible(PyObject* obj_ptr){
		// the second condition is important, for some reason otherwise there were attempted conversions of Body to list which failed afterwards.
		if(!PySequence_Check(obj_ptr) || !PyObject_HasAttrString(obj_ptr,"__len__")) return 0;
		return obj_ptr;
	}
	static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data){
		 void* storage=((boost::python::converter::rvalue_from_python_storage<std::vector<containedType> >*)(data))->storage.bytes;
		 new (storage) std::vector<containedType>();
		 std::vector<containedType>* v=(std::vector<containedType>*)(storage);
		 int l=PySequence_Size(obj_ptr); if(l<0) abort(); /*std::cerr<<"l="<<l<<"; "<<typeid(containedType).name()<<std::endl;*/ v->reserve(l); for(int i=0; i<l; i++) { v->push_back(boost::python::extract<containedType>(PySequence_GetItem(obj_ptr,i))); }
		 data->convertible=storage;
	}
};


struct custom_ptrMatchMaker_from_float{
	custom_ptrMatchMaker_from_float(){ boost::python::converter::registry::push_back(&convertible,&construct,boost::python::type_id<shared_ptr<MatchMaker> >()); }
	static void* convertible(PyObject* obj_ptr){ if(!PyNumber_Check(obj_ptr)) { cerr<<"Not convertible to MatchMaker"<<endl; return 0; } return obj_ptr; }
	static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data){
		void* storage=((boost::python::converter::rvalue_from_python_storage<shared_ptr<MatchMaker> >*)(data))->storage.bytes;
		new (storage) shared_ptr<MatchMaker>(new MatchMaker); // allocate the object at given address
		shared_ptr<MatchMaker>* mm=(shared_ptr<MatchMaker>*)(storage); // convert that address to our type
		(*mm)->algo="val"; (*mm)->val=PyFloat_AsDouble(obj_ptr); (*mm)->postLoad(**mm);
		data->convertible=storage;
	}
};



#ifdef YADE_MASK_ARBITRARY
struct custom_mask_to_long{
	static PyObject* convert(const mask_t& mask){
		return PyLong_FromString(const_cast<char*>(mask.to_string().c_str()),NULL,2);
	}
};
struct custom_mask_from_long{
	custom_mask_from_long(){
		 boost::python::converter::registry::push_back(&convertible,&construct,boost::python::type_id<mask_t>());
	}
	static void* convertible(PyObject* obj_ptr){
		return (PyLong_Check(obj_ptr) || PyInt_Check(obj_ptr))? obj_ptr : 0;
	}
	static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data){
		void* storage=((boost::python::converter::rvalue_from_python_storage<mask_t>*)(data))->storage.bytes;
		new (storage) mask_t; mask_t* mask=(mask_t*)storage;
		if (PyInt_Check(obj_ptr)) obj_ptr = PyLong_FromLong(PyInt_AsLong(obj_ptr));
		obj_ptr = _PyLong_Format(obj_ptr,2,0,0);
		std::string s(PyString_AsString(obj_ptr));
		//
		if (s.substr(0,2).compare("0b")==0) s = s.substr(2);
		if (s[s.length()-1]=='L') s = s.substr(0,s.length()-1);
		// TODO?
		*mask = mask_t(s);
		data->convertible=storage;
	}
};
#endif


BOOST_PYTHON_MODULE(_customConverters){

	custom_Se3r_from_seq(); boost::python::to_python_converter<Se3r,custom_se3_to_tuple>();

	custom_OpenMPAccumulator_from_float(); boost::python::to_python_converter<OpenMPAccumulator<Real>, custom_OpenMPAccumulator_to_float>(); 
	custom_OpenMPAccumulator_from_int(); boost::python::to_python_converter<OpenMPAccumulator<int>, custom_OpenMPAccumulator_to_int>(); 
	// todo: OpenMPAccumulator<int>

	custom_ptrMatchMaker_from_float();

	// StrArrayMap (typedef for std::map<std::string,numpy_boost>) → python dictionary
	//custom_StrArrayMap_to_dict();
	// register from-python converter and to-python converter

	boost::python::to_python_converter<std::vector<std::vector<std::string> >,custom_vvector_to_list<std::string> >();
	//boost::python::to_python_converter<std::list<shared_ptr<Functor> >, custom_list_to_list<shared_ptr<Functor> > >();
	//boost::python::to_python_converter<std::list<shared_ptr<Functor> >, custom_list_to_list<shared_ptr<Functor> > >();

#ifdef YADE_MASK_ARBITRARY
	custom_mask_from_long();
	boost::python::to_python_converter<mask_t,custom_mask_to_long>();
#endif

	// register 2-way conversion between c++ vector and python homogeneous sequence (list/tuple) of corresponding type
	#define VECTOR_SEQ_CONV(Type) custom_vector_from_seq<Type>();  boost::python::to_python_converter<std::vector<Type>, custom_vector_to_list<Type> >();
		VECTOR_SEQ_CONV(int);
		VECTOR_SEQ_CONV(bool);
		VECTOR_SEQ_CONV(Real);
		VECTOR_SEQ_CONV(Se3r);
		VECTOR_SEQ_CONV(Vector2r);
		VECTOR_SEQ_CONV(Vector2i);
		VECTOR_SEQ_CONV(Vector3r);
		VECTOR_SEQ_CONV(Vector3i);
		VECTOR_SEQ_CONV(Vector6r);
		VECTOR_SEQ_CONV(Vector6i);
		VECTOR_SEQ_CONV(Matrix3r);
		VECTOR_SEQ_CONV(Matrix6r);
		VECTOR_SEQ_CONV(std::string);
		VECTOR_SEQ_CONV(shared_ptr<Body>);
		VECTOR_SEQ_CONV(shared_ptr<Engine>);
		VECTOR_SEQ_CONV(shared_ptr<Material>);
		VECTOR_SEQ_CONV(shared_ptr<Serializable>);
		VECTOR_SEQ_CONV(shared_ptr<BoundFunctor>);
		VECTOR_SEQ_CONV(shared_ptr<IGeomFunctor>);
		VECTOR_SEQ_CONV(shared_ptr<IPhysFunctor>);
		VECTOR_SEQ_CONV(shared_ptr<LawFunctor>);
		VECTOR_SEQ_CONV(shared_ptr<IntrCallback>);
		#ifdef YADE_BODY_CALLBACK
			VECTOR_SEQ_CONV(shared_ptr<BodyCallback>);
		#endif
		VECTOR_SEQ_CONV(shared_ptr<SpherePack>);
		VECTOR_SEQ_CONV(shared_ptr<KinematicEngine>);
		#ifdef YADE_OPENGL
			VECTOR_SEQ_CONV(shared_ptr<GlBoundFunctor>);
			VECTOR_SEQ_CONV(shared_ptr<GlStateFunctor>);
			VECTOR_SEQ_CONV(shared_ptr<GlShapeFunctor>);
			VECTOR_SEQ_CONV(shared_ptr<GlIGeomFunctor>);
			VECTOR_SEQ_CONV(shared_ptr<GlIPhysFunctor>);
			VECTOR_SEQ_CONV(shared_ptr<GlExtraDrawer>);
		#endif
	#undef VECTOR_SEQ_CONV
}





