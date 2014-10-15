// © 2013 Jan Elias, http://www.fce.vutbr.cz/STM/elias.j/, elias.j@fce.vutbr.cz
// https://www.vutbr.cz/www_base/gigadisk.php?i=95194aa9a

#ifdef YADE_CGAL

#include"Polyhedra.hpp"
#include<pkg/common/Sphere.hpp>

//***************************************************************************
/*! Create Polyhedra (collision geometry) from colliding Polyhedras. */
class Ig2_Polyhedra_Polyhedra_PolyhedraGeom: public IGeomFunctor
{
	public:
		virtual ~Ig2_Polyhedra_Polyhedra_PolyhedraGeom(){};
		virtual bool go(const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c);
		virtual bool goReverse(	const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c) { return go(cm1,cm2,state2,state1,-shift2,force,c); }
		FUNCTOR2D(Polyhedra,Polyhedra);
		DEFINE_FUNCTOR_ORDER_2D(Polyhedra,Polyhedra);
		YADE_CLASS_BASE_DOC(Ig2_Polyhedra_Polyhedra_PolyhedraGeom,IGeomFunctor,"Create/update geometry of collision between 2 Polyhedras");	
		DECLARE_LOGGER;	
	private:
};
REGISTER_SERIALIZABLE(Ig2_Polyhedra_Polyhedra_PolyhedraGeom);

//***************************************************************************
/*! Create Polyhedra (collision geometry) from colliding Wall & Polyhedra. */
class Ig2_Wall_Polyhedra_PolyhedraGeom: public IGeomFunctor
{
	public:
		virtual ~Ig2_Wall_Polyhedra_PolyhedraGeom(){};
		virtual bool go(const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c);
		FUNCTOR2D(Wall,Polyhedra);
		DEFINE_FUNCTOR_ORDER_2D(Wall,Polyhedra);
		YADE_CLASS_BASE_DOC(Ig2_Wall_Polyhedra_PolyhedraGeom,IGeomFunctor,"Create/update geometry of collision between Wall and Polyhedra");	
		DECLARE_LOGGER;	
	private:		
};
REGISTER_SERIALIZABLE(Ig2_Wall_Polyhedra_PolyhedraGeom);

//***************************************************************************
/*! Create Polyhedra (collision geometry) from colliding Facet & Polyhedra. */
class Ig2_Facet_Polyhedra_PolyhedraGeom: public IGeomFunctor
{
	public:
		virtual ~Ig2_Facet_Polyhedra_PolyhedraGeom(){};
		virtual bool go(const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c);
		FUNCTOR2D(Facet,Polyhedra);
		DEFINE_FUNCTOR_ORDER_2D(Facet,Polyhedra);
		YADE_CLASS_BASE_DOC(Ig2_Facet_Polyhedra_PolyhedraGeom,IGeomFunctor,"Create/update geometry of collision between Facet and Polyhedra");	
		DECLARE_LOGGER;	
	private:		
};
REGISTER_SERIALIZABLE(Ig2_Facet_Polyhedra_PolyhedraGeom);

//***************************************************************************
/*! Create Polyhedra (collision geometry) from colliding Sphere & Polyhedra. */
/*
class Ig2_Sphere_Polyhedra_PolyhedraGeom: public IGeomFunctor
{
	public:
		virtual ~Ig2_Sphere_Polyhedra_PolyhedraGeom(){};
		virtual bool go(const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c);
		FUNCTOR2D(Sphere,Polyhedra);
		DEFINE_FUNCTOR_ORDER_2D(Sphere,Polyhedra);
		YADE_CLASS_BASE_DOC(Ig2_Sphere_Polyhedra_PolyhedraGeom,IGeomFunctor,"Create/update geometry of collision between Sphere and Polyhedra");	
		DECLARE_LOGGER;	
	private:		
};
REGISTER_SERIALIZABLE(Ig2_Sphere_Polyhedra_PolyhedraGeom);
*/

#endif // YADE_CGAL
