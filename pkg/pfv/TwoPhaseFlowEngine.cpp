 
/*************************************************************************
*  Copyright (C) 2014 by Bruno Chareyre <bruno.chareyre@hmg.inpg.fr>     *
*  Copyright (C) 2013 by T. Sweijen (T.sweijen@uu.nl)                    *
*  Copyright (C) 2012 by Chao Yuan <chao.yuan@3sr-grenoble.fr>           *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

// This is an example of how to derive a new FlowEngine with additional data and possibly completely new behaviour.
// Every functions of the base engine can be overloaded, and new functions can be added

//keep this #ifdef as long as you don't really want to realize a final version publicly, it will save compilation time for everyone else
//when you want it compiled, you can pass -DDFNFLOW to cmake, or just uncomment the following line
#define TWOPHASEFLOW
#ifdef TWOPHASEFLOW

#include "FlowEngine_TwoPhaseFlowEngineT.hpp"

/// We can add data to the Info types by inheritance
class TwoPhaseCellInfo : public FlowCellInfo_TwoPhaseFlowEngineT
{
	public:
  	bool isWRes;
	bool isNWRes;
	bool isTrapW;
	bool isTrapNW;
	double saturation;
	bool isImbibition;
	double trapCapP;//for calculating the pressure of trapped phase, cell->info().p() = pressureNW- trapCapP. OR cell->info().p() = pressureW + trapCapP
	std::vector<double> poreThroatRadius;
	double poreBodyRadius;
	double poreBodyVolume;
	TwoPhaseCellInfo (void)
	{
		isWRes = true; isNWRes = false; isTrapW = false; isTrapNW = false;
		saturation = 1.0;
		isImbibition = false;
		trapCapP = 0;
		poreThroatRadius.resize(4, 0);
		poreBodyRadius = 0;
		poreBodyVolume = 0;
	}
	
};

class TwoPhaseVertexInfo : public FlowVertexInfo_TwoPhaseFlowEngineT {
	public:
	//same here if needed
};

typedef TemplateFlowEngine_TwoPhaseFlowEngineT<TwoPhaseCellInfo,TwoPhaseVertexInfo> TwoPhaseFlowEngineT;
REGISTER_SERIALIZABLE(TwoPhaseFlowEngineT);
YADE_PLUGIN((TwoPhaseFlowEngineT));

class TwoPhaseFlowEngine : public TwoPhaseFlowEngineT
{
		double totalCellVolume;
	public :
	//We can overload every functions of the base engine to make it behave differently
	//if we overload action() like this, this engine is doing nothing in a standard timestep, it can still have useful functions
	virtual void action() {};
	
	//If a new function is specific to the derived engine, put it here, else go to the base TemplateFlowEngine
	//if it is useful for everyone
	void fancyFunction(Real what);
	void computePoreBodyRadius();
	YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(TwoPhaseFlowEngine,TwoPhaseFlowEngineT,"documentation here",
	((double,surfaceTension,0.0728,,"Water Surface Tension in contact with air at 20 Degrees Celsius is: 0.0728(N/m)"))
	((bool,initialWetting,true,,"Initial wetting saturated (=true) or non-wetting saturated (=false)"))
	((bool, isPhaseTrapped,true,,"If True, both phases can be entrapped by the other, which would correspond to snap-off. If false, both phases are always connected to their reservoirs, thus no snap-off."))
	
	,/*TwoPhaseFlowEngineT()*/,
	,
	.def("fancyFunction",&TwoPhaseFlowEngine::fancyFunction,(boost::python::arg("what")=0),"test function")
	)
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(TwoPhaseFlowEngine);
YADE_PLUGIN((TwoPhaseFlowEngine));

void TwoPhaseFlowEngine::fancyFunction(Real what) {std::cerr<<"yes, I'm a new function"<<std::endl;}


void TwoPhaseFlowEngine::computePoreBodyRadius()
{
  
    // This routine finds the radius of the inscribed sphere within each pore-body
    // Following Mackay et al., 1972. 
      double d01 = 0.0, d02 = 0.0, d03 = 0.0, d12 = 0.0, d13 = 0.0, d23 = 0.0, Rin = 0.0, r0 = 0.0, r1 = 0.0, r2 =0.0, r3 = 0.0;
      bool check = false;
      unsigned int i = 0;
      double dR=0.0, tempR = 0.0;
      bool initialSign = false; //False = negative, true is positive
      bool first = true;
      
      Eigen::MatrixXd M(6,6);


      
      FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles){
	//Distance between multiple particles, can be done more efficient
	d01 = d02 = d03 = d12 = d13 = d23 = r0 = r1 = r2= r3 = 0.0;
	
	d01 = pow((cell->vertex(0)->point().x()-cell->vertex(1)->point().x()),2)+
	pow((cell->vertex(0)->point().y()-cell->vertex(1)->point().y()),2)+
	pow((cell->vertex(0)->point().z()-cell->vertex(1)->point().z()),2);
	
	d02 = pow((cell->vertex(0)->point().x()-cell->vertex(2)->point().x()),2)+
	pow((cell->vertex(0)->point().y()-cell->vertex(2)->point().y()),2)+
	pow((cell->vertex(0)->point().z()-cell->vertex(2)->point().z()),2);
	
	d03 = pow((cell->vertex(0)->point().x()-cell->vertex(3)->point().x()),2)+
	pow((cell->vertex(0)->point().y()-cell->vertex(3)->point().y()),2)+
	pow((cell->vertex(0)->point().z()-cell->vertex(3)->point().z()),2);
	
	d12 =pow((cell->vertex(1)->point().x()-cell->vertex(2)->point().x()),2)+
	pow((cell->vertex(1)->point().y()-cell->vertex(2)->point().y()),2)+
	pow((cell->vertex(1)->point().z()-cell->vertex(2)->point().z()),2);
	
	d13 = pow((cell->vertex(1)->point().x()-cell->vertex(3)->point().x()),2)+
	pow((cell->vertex(1)->point().y()-cell->vertex(3)->point().y()),2)+
	pow((cell->vertex(1)->point().z()-cell->vertex(3)->point().z()),2);
	
	d23 = pow((cell->vertex(2)->point().x()-cell->vertex(3)->point().x()),2)+
	pow((cell->vertex(2)->point().y()-cell->vertex(3)->point().y()),2)+
	pow((cell->vertex(2)->point().z()-cell->vertex(3)->point().z()),2);
	

	//Radii of the particles
	r0 = sqrt(cell -> vertex(0) -> point().weight());
	r1 = sqrt(cell -> vertex(1) -> point().weight());
	r2 = sqrt(cell -> vertex(2) -> point().weight());
	r3 = sqrt(cell -> vertex(3) -> point().weight());
	
	
	//Fill coefficient matrix
	M(0,0) = 0.0;
	M(1,0) = d01;
	M(2,0) = d02;
	M(3,0) = d03;
	M(4,0) = pow((r0+Rin),2);
	M(5,0) = 1.0;
	
	M(0,1) = d01;
	M(1,1) = 0.0;
	M(2,1) = d12;
	M(3,1) = d13;
	M(4,1) = pow((r1+Rin),2);
	M(5,1) = 1.0;
	
	M(0,2) = d02; 
	M(1,2) = d12;
	M(2,2) = 0.0;
	M(3,2) = d23;
	M(4,2) = pow((r2+Rin),2);
	M(5,2) = 1.0;
	
	M(0,3) = d03;
	M(1,3) = d13;
	M(2,3) = d23;
	M(3,3) = 0.0;
	M(4,3) = pow((r3+Rin),2);
	M(5,3) = 1.0;
	
	M(0,4) = pow((r0+Rin),2);
	M(1,4) = pow((r1+Rin),2);
	M(2,4) = pow((r2+Rin),2);
	M(3,4) = pow((r3+Rin),2);
	M(4,4) = 0.0;
	M(5,4) = 1.0;
	
	M(0,5) = 1.0;
	M(1,5) = 1.0;
	M(2,5) = 1.0;
	M(3,5) = 1.0;
	M(4,5) = 1.0;
	M(5,5) = 0.0;
	

	i = 0;
	check = false;
	dR = Rin = 0.0 + (min(r0,min(r1,min(r2,r3))) / 50.0); //Estimate an initial dR
	first = true;
	//Iterate untill check = true, such that an accurate answer as been found
	while (check == false){
	  i = i + 1;
	  tempR = Rin;
	  Rin = Rin + dR;
	
	  M(4,0) = pow((r0+Rin),2);
	  M(4,1) = pow((r1+Rin),2);
	  M(4,2) = pow((r2+Rin),2);
	  M(4,3) = pow((r3+Rin),2);
	  M(0,4) = pow((r0+Rin),2);
	  M(1,4) = pow((r1+Rin),2);
	  M(2,4) = pow((r2+Rin),2);
	  M(3,4) = pow((r3+Rin),2);
	
	
	  if (first){
	    first = false;
	    if(M.determinant() < 0.0){initialSign = false;} //Initial D is negative
	    if(M.determinant() > 0.0){initialSign = true;} // Initial D is positive
	  }
	
	  if(std::abs(M.determinant()) < 1E-30){check = true;}
	
	
	  if((initialSign==true) && (check ==false)){
	    if(M.determinant() < 0.0){
	      Rin = Rin -dR;
	      dR = dR / 2.0; 
	    }	  
	  }
	
	  if((initialSign==false) && (check ==false)){
	    if(M.determinant() > 0.0){
	      Rin = Rin -dR;
	      dR = dR / 2.0;  
	    }	  
	  }
	
	  cout << endl << i << " "<<Rin << " "<< dR << " "<< M.determinant();
	  if(i > 4000){
	    cout << endl << "error, finding solution takes too long cell:" << cell->info().id;
	    check = true;
	  }
	  if ( std::abs(tempR - Rin)/Rin < 0.001){check = true;}
	
      }
    cell -> info().poreBodyRadius = Rin;
   }
}




#endif //TwoPhaseFLOW
 
