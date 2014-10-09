/*************************************************************************
*  Copyright (C) 2013 by T. Sweijen (T.sweijen@uu.nl)                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#ifdef YADE_CGAL
#ifdef FLOW_ENGINE

#define SOLUTE_FLOW
#ifdef SOLUTE_FLOW

#include "FlowEngine_SoluteFlowEngineT.hpp"

#include <Eigen/Sparse>
#include <Eigen/Dense>
// #include "FlowEngine.hpp.in"

class SoluteCellInfo : public FlowCellInfo_SoluteFlowEngineT
{	
	public:
 	Real WSaturation2;
	Real InsideSphereRadius2;
	bool imbibition2;
	bool networkConnectivity2;
	bool NnetworkConnectivity2;
	bool error2;
	std::vector<double> inscribedCircleRadius2;
	
	
// 	liquidbr.resize(4, 0);inline std::vector<double>& Rh (void) {return rayHydr;}
 	SoluteCellInfo (void) : FlowCellInfo_SoluteFlowEngineT() {WSaturation2=0.0;InsideSphereRadius2=0.0;inscribedCircleRadius2.resize(4,0);}
 	inline Real& saturationTEMP (void) {return WSaturation2;}
 	inline const Real& saturationTEMP (void) const {return WSaturation2;}
	inline std::vector<double>& inscribedCircleRadius (void) {return inscribedCircleRadius2;}
 	
 	
  	inline Real& InsideSphereRadius (void) {return InsideSphereRadius2;}
  	inline const Real& InsideSphereRadius (void) const {return InsideSphereRadius2;}
  	inline bool& imbibition (void) {return imbibition2;} //If imbibition is true, an interface is present within pore body
  	inline const bool& imbibition (void) const {return imbibition2;}
  	inline bool& Wnetwork (void) {return networkConnectivity2;} //If imbibition is true, an interface is present within pore body
  	inline const bool& Wnetwork (void) const {return networkConnectivity2;}
  	inline bool& error (void) {return error2;} //If imbibition is true, an interface is present within pore body
  	inline const bool& error (void) const {return error2;}
  	inline bool& NWnetwork (void) {return NnetworkConnectivity2;} //If imbibition is true, an interface is present within pore body
  	inline const bool& NWnetwork (void) const {return NnetworkConnectivity2;}
 	inline void getInfo (const SoluteCellInfo& otherCellInfo) {FlowCellInfo_SoluteFlowEngineT::getInfo(otherCellInfo); saturationTEMP()=otherCellInfo.saturationTEMP();}
};

typedef TemplateFlowEngine_SoluteFlowEngineT<SoluteCellInfo,FlowVertexInfo_SoluteFlowEngineT> SoluteFlowEngineT;
REGISTER_SERIALIZABLE(SoluteFlowEngineT);
YADE_PLUGIN((SoluteFlowEngineT));

class SoluteFlowEngine : public SoluteFlowEngineT
{
	public :
		void getInsideEqSphere();
		void drainageFunction();
		void TPFaction();
// 		void LiquidBridge();
// 		double liquidbridgeCalc(double R1, double R2, double Pc);
		void copySaturation();
		void copySaturationTEMP();
		void inscribedCircle();
		void staticIterator();
		void InitializeScenario();
		void imbibitionFunction();
		void InsertBoundaryConditions();
		void imbibitionPoreFunction();
		void resetNetwork();
		void checkConnectivityNetwork();
		void PorethroatSaturation(double Pc, int ID,unsigned int facet);
		void saveVtk(const char* folder);
		double getSaturation(unsigned int ID){return solver->T[solver->currentTes].cellHandles[ID]->info().saturation();}
		double getInscribedRadius(unsigned int ID){return solver->T[solver->currentTes].cellHandles[ID]->info().InsideSphereRadius();}

		double getBulkSaturation();

		///Elaborate the description as you wish
		YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(SoluteFlowEngine,SoluteFlowEngineT,"A variant of :yref:`FlowEngine` with solute transport).",
		///No additional variable yet, else input here
		((double,initialWettingSaturation,1.0,,"Initial water saturation inside the model"))
		((double,InterfacialTension,72.0,,"Interfacial tension between two phases in [dynes/cm]"))
		((double,WettingAnlge,0.0,,"Wetting angle in degrees."))
		((double,YoungsLaplaceConstante,0.0,,"2*interfacial tension* cos(wetting angle)"))
		((unsigned int,BCIDWettingPhase,2,,"Boundary ID of the wetting Phase (default =3)"))
		((double,BCPressureWettingPhase,0.0,,"Pressure at the Boundary of the Wetting Phase"))
		((unsigned int,BCIDNonWettingPhase,3,,"Boundary ID of the Non-wetting Phase (default =2)"))
		((double,BCPressureNonWettingPhase,100.0,,"Pressure at the Boundary of the Non-Wetting Phase"))
		((double,CappilaryPressureOverall,0.0,,"Difference in phase pressures as found at the boundaries"))
		((bool,first,true,,"First Calculation? "))
		((bool,active,true,,"Is something still happening?"))
		((bool,makeMovie,true,,"Make a movie of the different steps to reach equilibrium?"))
		((bool,drainage,true,,"Do we have drainage or nor?"))
		((bool,activeiterator,true,,"Do we have drainage or nor?"))
		((bool,cornerflow,false,,"Do we have corner flow during imbibition or not?"))

		,,,
		.def("TPFaction",&SoluteFlowEngine::TPFaction,"Activate Two Phase Flow Engine")
 		.def("PorethroatSaturation",&SoluteFlowEngine::PorethroatSaturation,(boost::python::arg("Pc"),boost::python::arg("ID"),boost::python::arg("facet")),"Calculate swelling of liquid bridge")
		.def("getBulkSaturation",&SoluteFlowEngine::getBulkSaturation,"Activate Two Phase Flow Engine")
		.def("staticIterator",&SoluteFlowEngine::staticIterator,"Iterator to find eq. in two-phase flow")
		.def("getSaturation",&SoluteFlowEngine::getSaturation,(boost::python::arg("ID")),"Get Saturation")
		.def("getInscribedRadius",&SoluteFlowEngine::getInscribedRadius,(boost::python::arg("ID")),"Get Saturation")
		.def("getInsideEqSphere",&SoluteFlowEngine::getInsideEqSphere,"Get equivalent pore volume")
		.def("InitializeScenario",&SoluteFlowEngine::InitializeScenario,"Reset Initial conditions")
		.def("InsertBoundaryConditions",&SoluteFlowEngine::InsertBoundaryConditions,"Insert the boundary conditions into vectors.")
		.def("checkConnectivityNetwork",&SoluteFlowEngine::checkConnectivityNetwork,"Reset Initial conditions")
		.def("inscribedCircle",&SoluteFlowEngine::inscribedCircle,"calculate the radius of an inscribed circle")

	
		)
};
REGISTER_SERIALIZABLE(SoluteFlowEngine);


// PeriodicFlowEngine::~PeriodicFlowEngine(){}


void SoluteFlowEngine::PorethroatSaturation(double Pc, int ID, unsigned int facet)
{
         // NOTE: Instead of passing on ID, I could pass on a CellHandle& cell, but that is for internal within C++
        FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles){
	  if (cell->info().id == ID){

        double R0=0.0, R1=0.0, R2=0.0, R3=0.0;

     
	  Eigen::MatrixXf M(2,3), C(2,3);
	 Eigen::VectorXf X(3), Z(3);
	 Eigen::VectorXf B(2), D(2);
 	  M(0,0) = 2*(cell->vertex(facetVertices[facet][0])->point().x() -  cell->vertex(facetVertices[facet][1])->point().x());
 	  M(1,0) = 2*(cell->vertex(facetVertices[facet][0])->point().x() -  cell->vertex(facetVertices[facet][2])->point().x());
	  M(0,1) = 2*(cell->vertex(facetVertices[facet][0])->point().y() -  cell->vertex(facetVertices[facet][1])->point().y());
	  M(1,1) = 2*(cell->vertex(facetVertices[facet][0])->point().y() -  cell->vertex(facetVertices[facet][2])->point().y());
	  M(0,2) = 2*(cell->vertex(facetVertices[facet][0])->point().z() -  cell->vertex(facetVertices[facet][1])->point().z());
	  M(1,2) = 2*(cell->vertex(facetVertices[facet][0])->point().z() -  cell->vertex(facetVertices[facet][2])->point().z());


	  
	  B(0) = (pow(cell->vertex(facetVertices[facet][0])->point().x(),2)+pow(cell->vertex(facetVertices[facet][0])->point().y(),2)
	      +pow(cell->vertex(facetVertices[facet][0])->point().z(),2) - cell->vertex(facetVertices[facet][0])->point().weight())
	      -(pow(cell->vertex(facetVertices[facet][1])->point().x(),2)+pow(cell->vertex(facetVertices[facet][1])->point().y(),2)
	       +pow(cell->vertex(facetVertices[facet][1])->point().z(),2) - cell->vertex(facetVertices[facet][1])->point().weight())
	      -2*Pc*(sqrt(cell->vertex(facetVertices[facet][0])->point().weight())-sqrt(cell->vertex(facetVertices[facet][1])->point().weight()));
	  
	  B(1) = (pow(cell->vertex(facetVertices[facet][0])->point().x(),2)+pow(cell->vertex(facetVertices[facet][0])->point().y(),2)	  
	      +pow(cell->vertex(facetVertices[facet][0])->point().z(),2) - cell->vertex(facetVertices[facet][0])->point().weight())	  
	      -(pow(cell->vertex(facetVertices[facet][2])->point().x(),2)+pow(cell->vertex(facetVertices[facet][2])->point().y(),2)	  
	       +pow(cell->vertex(facetVertices[facet][2])->point().z(),2) - cell->vertex(facetVertices[facet][2])->point().weight())	  
	      -2*Pc*(sqrt(cell->vertex(facetVertices[facet][0])->point().weight())-sqrt(cell->vertex(facetVertices[facet][1])->point().weight()));	  
	
	    X = M.fullPivLu().solve(B);  


	    cout << endl <<B[0]<< " " << B[1];
	      cout << endl << X[0] << " "<< X[1] << " "<< X[2];
	  
       cout << "hallo!";
       
       R0= sqrt(cell->vertex(0)->point().weight());
       R1 = sqrt(cell->vertex(1)->point().weight());
       R2 = sqrt(cell->vertex(2)->point().weight());
       R3 = sqrt(cell->vertex(3)->point().weight());
       
       
       
	  }
	}
  
}








void SoluteFlowEngine::TPFaction()
{
  
  // Initialize calculations, add Boundary conditions and calculate different constants
  if(first){
    InitializeScenario();
    InsertBoundaryConditions();
    first = false;
    cout << endl << "Calculate radius of inscribed sphere";
    getInsideEqSphere();
     //inscribedCircle();
    cout << endl << "Calculate radius of inscribed sphere is done!";
    }
  
  activeiterator = true;
  if(cornerflow == true){resetNetwork();}
  // While loop for calculating equilibrium conditions for a set Pn-Pw, only solved for continuous flow. 
   while(activeiterator){
  if(cornerflow == false){checkConnectivityNetwork();}
  copySaturation();
  activeiterator = false;
  InsertBoundaryConditions();
  staticIterator();
  copySaturationTEMP();
  if(makeMovie){solver->saveVtk("./VTK");}  
  }

  
  
}

void SoluteFlowEngine::staticIterator()
{
    //This is an intermediate state function, merely do distinquish between drainage and imbibition
	if(drainage){drainageFunction();}
	if(drainage == false){
	  imbibitionFunction();
	  imbibitionPoreFunction();
	}
} 











void SoluteFlowEngine::inscribedCircle()
{
  //This function calculates the inscribed circle of a pore throat facet, for non-touching spheres, al though it yields then 
  
  //same results as the effectiveradius in solver.
  double d01 = 0.0, d02 = 0.0, d12 = 0.0,d03 = 0.0, d13 = 0.0, d23 = 0.0;
  double R0 = 0.0, R1 = 0.0, R2 = 0.0, R3 = 0.0, Rm0 = 0.0, Rm1 = 0.0, Rm2 = 0.0, Rm3 = 0.0;
  double cosalpha0 = 0.0, cosalpha1 = 0.0, cosalpha2 = 0.0, cosalpha3 = 0.0;
  double cosalpha10 = 0.0, cosalpha11 = 0.0,cosalpha12 = 0.0, cosalpha13 = 0.0;
  double cosalpha20 = 0.0, cosalpha21 = 0.0,cosalpha22 = 0.0, cosalpha23 = 0.0;
  FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles){
  R0= sqrt(cell->vertex(0)->point().weight());
  R1 = sqrt(cell->vertex(1)->point().weight());
  R2 = sqrt(cell->vertex(2)->point().weight());
  R3 = sqrt(cell->vertex(3)->point().weight());
  
  	d01 = sqrt(pow((cell->vertex(0)->point().x()-cell->vertex(1)->point().x()),2)+
	pow((cell->vertex(0)->point().y()-cell->vertex(1)->point().y()),2)+
	pow((cell->vertex(0)->point().z()-cell->vertex(1)->point().z()),2));
	
	d02 = sqrt(pow((cell->vertex(0)->point().x()-cell->vertex(2)->point().x()),2)+
	pow((cell->vertex(0)->point().y()-cell->vertex(2)->point().y()),2)+
	pow((cell->vertex(0)->point().z()-cell->vertex(2)->point().z()),2));
	
	d03 = sqrt(pow((cell->vertex(0)->point().x()-cell->vertex(3)->point().x()),2)+
	pow((cell->vertex(0)->point().y()-cell->vertex(3)->point().y()),2)+
	pow((cell->vertex(0)->point().z()-cell->vertex(3)->point().z()),2));
	
	d12 =sqrt(pow((cell->vertex(1)->point().x()-cell->vertex(2)->point().x()),2)+
	pow((cell->vertex(1)->point().y()-cell->vertex(2)->point().y()),2)+
	pow((cell->vertex(1)->point().z()-cell->vertex(2)->point().z()),2));
	
	d13 = sqrt(pow((cell->vertex(1)->point().x()-cell->vertex(3)->point().x()),2)+
	pow((cell->vertex(1)->point().y()-cell->vertex(3)->point().y()),2)+
	pow((cell->vertex(1)->point().z()-cell->vertex(3)->point().z()),2));
	
	d23 = sqrt(pow((cell->vertex(2)->point().x()-cell->vertex(3)->point().x()),2)+
	pow((cell->vertex(2)->point().y()-cell->vertex(3)->point().y()),2)+
	pow((cell->vertex(2)->point().z()-cell->vertex(3)->point().z()),2));
	
	cosalpha0 = (pow(d23,2)-pow(d12,2)-pow(d13,2)) / (-2.0*d12*d13);
	cosalpha1 = (pow(d23,2)-pow(d02,2)-pow(d03,2)) / (-2.0*d02*d03);
	cosalpha2 = (pow(d03,2)-pow(d01,2)-pow(d13,2)) / (-2.0*d01*d13);
	cosalpha3 = (pow(d02,2)-pow(d12,2)-pow(d01,2)) / (-2.0*d12*d01);




	

	 Rm0 = Rm1 = Rm2 = Rm3 = R1;
	  for(unsigned int i=0; i<2000; i++){
  
	cosalpha10 = (pow((R2+Rm0),2)-pow((R1+Rm0),2)-pow(d12,2))/(-2*d12*(R1+Rm0));
	cosalpha11 = (pow((R2+Rm1),2)-pow((R0+Rm1),2)-pow(d02,2))/(-2*d12*(R0+Rm1));
	cosalpha12 = (pow((R0+Rm2),2)-pow((R1+Rm2),2)-pow(d01,2))/(-2*d01*(R1+Rm2));
	cosalpha13 = (pow((R2+Rm3),2)-pow((R1+Rm3),2)-pow(d12,2))/(-2*d12*(R1+Rm3));
	
	cosalpha20 = (pow((R3+Rm0),2)-pow((R1+Rm0),2)-pow(d13,2))/(-2*d13*(R1+Rm0));
	cosalpha21 = (pow((R3+Rm1),2)-pow((R0+Rm1),2)-pow(d03,2))/(-2*d03*(R0+Rm1));
	cosalpha22 = (pow((R3+Rm2),2)-pow((R1+Rm2),2)-pow(d13,2))/(-2*d13*(R1+Rm2));
	cosalpha23 = (pow((R0+Rm3),2)-pow((R1+Rm3),2)-pow(d01,2))/(-2*d01*(R1+Rm3));
	
	 if ( (acos(cosalpha0) - acos(cosalpha10) - acos(cosalpha20)) < 0.0){Rm0 = Rm0*0.99;}
	 if ( (acos(cosalpha0) - acos(cosalpha10) - acos(cosalpha20)) == 0.0){Rm0 = Rm0;}
	 if ( (acos(cosalpha0) - acos(cosalpha10) - acos(cosalpha20)) > 0.0){Rm0 = Rm0*1.01;}
	    
	 if ( (acos(cosalpha1) - acos(cosalpha11) - acos(cosalpha21)) < 0.0){Rm1 = Rm1*0.99;}	    
	 if ( (acos(cosalpha1) - acos(cosalpha11) - acos(cosalpha21)) == 0.0){Rm1 = Rm1;}	    
	 if ( (acos(cosalpha1) - acos(cosalpha11) - acos(cosalpha21)) > 0.0){Rm1 = Rm1*1.01;}	 
	 
	 if ( (acos(cosalpha2) - acos(cosalpha12) - acos(cosalpha22)) < 0.0){Rm2 = Rm2*0.99;}	 
	 if ( (acos(cosalpha2) - acos(cosalpha12) - acos(cosalpha22)) == 0.0){Rm2 = Rm2;}	 
	 if ( (acos(cosalpha2) - acos(cosalpha12) - acos(cosalpha22)) > 0.0){Rm2 = Rm2*1.01;}	 
	 
	 if ( (acos(cosalpha3) - acos(cosalpha13) - acos(cosalpha23)) < 0.0){Rm3 = Rm3*0.99;}	 
	 if ( (acos(cosalpha3) - acos(cosalpha13) - acos(cosalpha23)) == 0.0){Rm3 = Rm3;}	 
	 if ( (acos(cosalpha3) - acos(cosalpha13) - acos(cosalpha23)) > 0.0){Rm3 = Rm3*1.01;}	 
 
	  }
	  //cout << endl << Rm0 << " "<< Rm1 << " "<< Rm2 << " "<<Rm3;
	  cell->info().inscribedCircleRadius()[0] = Rm0;
	  cell->info().inscribedCircleRadius()[1] = Rm1;
	  cell->info().inscribedCircleRadius()[2] = Rm2;
	  cell->info().inscribedCircleRadius()[3] = Rm3;
	  
	  if(Rm0 > 1000.0*R1){ cell->info().inscribedCircleRadius()[0] = 0.0;}
	    if(Rm1 > 1000.0*R1){ cell->info().inscribedCircleRadius()[1] = 0.0;}
	    if(Rm2 > 1000.0*R1){ cell->info().inscribedCircleRadius()[2] = 0.0;}
	    if(Rm3 > 1000.0*R1){ cell->info().inscribedCircleRadius()[3] = 0.0;}
	  cout << endl << Rm0 << " "<< solver->computeEffectiveRadius(cell,0);
  
  }
  
}



void SoluteFlowEngine::resetNetwork()
{
    FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles)
 	{
	 cell->info().Wnetwork() = true;
	 cell->info().NWnetwork() = true;
 	}
  
}

void SoluteFlowEngine::InitializeScenario()
{
  FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles)
 	{
 	 cell->info().saturationTEMP() = initialWettingSaturation;
	 cell->info().saturation() = initialWettingSaturation;
	 cell->info().imbibition() = false;
	 cell->info().Wnetwork() = true;
	 cell->info().NWnetwork() = true;
	 cell->info().error() = false;
 	}

}


void SoluteFlowEngine::InsertBoundaryConditions()
{
    //Insert boundary condition saturations, pressure of the phase is assumed to be constant within on connected phase
    FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles)
    {
 		for (unsigned int ngb=0;ngb<4;ngb++){ 
 			if (cell->vertex(ngb)->info().id() == BCIDWettingPhase){
			  cell->info().saturation() = 1.0;
			  cell->info().saturationTEMP() = 1.0;
			  }

			if (cell->vertex(ngb)->info().id() == BCIDNonWettingPhase){
			  cell->info().saturation() = 0.0;
			  cell->info().saturationTEMP() = 0.0;
			}
 		}
    }
  
}
























void SoluteFlowEngine::drainageFunction()
{
   //Check pore-throat and inscribed circle for drainage.
    double Rcp = 0.0, throatRadius = 0.0;
    double CapillaryPressure = 0.0;
	CapillaryPressure = BCPressureNonWettingPhase - BCPressureWettingPhase;
	CappilaryPressureOverall = CapillaryPressure;
	Rcp = ((72.0 / 1000.0)) / CapillaryPressure;
	
    FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles)
    {
     if(cell->info().saturation() != 1.0){ 
     // Cell i is non-wetting saturated wereas the neighbor is water saturated 
      
	for (unsigned int ngb=0;ngb<4;ngb++){
	  if((cell->neighbor(ngb)->info().saturation() == 1.0)&&(cell->neighbor(ngb)->info().Wnetwork() == true)){
       
	  throatRadius = solver->computeEffectiveRadius(cell,ngb);//cell->info().inscribedCircleRadius()[ngb] //solver->computeEffectiveRadius(cell,ngb);
	  
	  //std::abs(solver->computeEffectiveRadius(cell, ngb));
	    if (throatRadius > Rcp){
		cell -> neighbor(ngb)->info().saturation() = 0.0; //NOTE(Thomas): include residual liquid bridge
     }
     }
     }
     }
  
   }

}
void SoluteFlowEngine::imbibitionPoreFunction()
{
       //Check hanging interface within a pore body
  	double Rcp = 0.0;
	double CapillaryPressure = 0.0;
	     CapillaryPressure = BCPressureNonWettingPhase - BCPressureWettingPhase;
	     CappilaryPressureOverall = CapillaryPressure;
	     Rcp = ((72.0 / 1000.0)) / CapillaryPressure;
	     
	FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles){
	  if(cell->info().imbibition() == true){
	    if(Rcp > cell->info().InsideSphereRadius()){
	      cell->info().imbibition() = false;
	      cell ->info().saturation() = 1.0; //NOTE(Thomas): include residual liquid bridge
	      
	   
	    
	    }
	  }
	}
	   
  
}



void SoluteFlowEngine::imbibitionFunction()
{
        //Reverse of drainage, and the option to assign hanging interfaces within pore-bodies.
	double Rcp = 0.0, throatRadius = 0.0;
	double CapillaryPressure = 0.0;
	  CapillaryPressure = BCPressureNonWettingPhase - BCPressureWettingPhase;
	  CappilaryPressureOverall = CapillaryPressure;
	  Rcp = ((72.0 / 1000.0)) / CapillaryPressure;
  
	  
	FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles)
	{
	    
	    if((cell->info().saturation() == 1.0)&&(cell->info().imbibition()==false)){ 
	      // Cell i is wetting saturated wereas the neighbor is non-wetting saturated 
	      for (unsigned int ngb=0;ngb<4;ngb++){
		if((cell->neighbor(ngb)->info().saturation() != 1.0) &&(cell->neighbor(ngb)->info().NWnetwork() == true)&&(cell->info().Wnetwork() == true)){
		  throatRadius = solver->computeEffectiveRadius(cell,ngb);//solver->computeEffectiveRadius(cell,ngb); //cell->info().inscribedCircleRadius()[ngb];
		    if (throatRadius <= Rcp){
		      cell -> neighbor(ngb)->info().imbibition() = true; // A new interface in ngb pore has been created
   
     }
     }

     }
  
   }
	  
	  
	  
  
}
}




void SoluteFlowEngine::checkConnectivityNetwork()
{
  bool networkactive;
  int temp;
   FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles)
   {
     cell->info().NWnetwork() = false;
     cell->info().Wnetwork() = false;
     
    //NOTE(Thomas): Check boundary conditions, may be implementen in insertBoundaryConditions
    for(unsigned vertex=0;vertex<4;vertex++){
      //NW connectivity
     if(cell->vertex(vertex)->info().id() ==  BCIDNonWettingPhase)
     {
       cell -> info().NWnetwork() = true; 
     }
     //W connectivity
     if(cell->vertex(vertex)->info().id() ==  BCIDWettingPhase)
     {
       cell -> info().Wnetwork() = true;
       
     }
    }
   }
 
    networkactive = true;
    temp = 0;
   while(networkactive){
     networkactive=false;
    FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles){
      if(cell -> info().NWnetwork() == true){
	for(unsigned ngb=0;ngb<4;ngb++){
	  if(cell->neighbor(ngb)->info().saturation() == 0.0){
	    networkactive = true;
	    cell->neighbor(ngb)->info().NWnetwork() = true;
	  }
	}
      }
      
      if(cell -> info().Wnetwork() == true){
	for(unsigned ngb=0;ngb<4;ngb++){
	  if(cell->neighbor(ngb)->info().saturation() == 1.0){
	    networkactive=true;
	    cell->neighbor(ngb)->info().Wnetwork() = true;
	  }
	}
      }
    }
    temp = temp+1;
    if(temp > 2000){networkactive = false;} // put upperlimit for while loop
    cout << endl << temp;
   }
   
   
    int Wcount = 0, NWcount = 0;
    FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles)
    {
      if(cell -> info().Wnetwork() == true){Wcount=Wcount+1;}
      if(cell -> info().NWnetwork() == true){NWcount=NWcount+1;}
    }
 cout << "NW " << NWcount << " W "<< Wcount;
}





double SoluteFlowEngine::getBulkSaturation()
{
     // This function calculates the bulk saturation of the wetting phase (Sw=1-Snw), using the summation of volumes. 

     double Vt = 0.0;
     double Vw = 0.0;
     bool boundary = false;
     FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles)
     {
       boundary = false;
       //Do not use the boundary cells.
       for (unsigned int ngb=0;ngb<4;ngb++){ 
	 if(cell->vertex(ngb)->info().id() < 6){boundary = true;}
	 
//  			if (cell->vertex(ngb)->info().id() == BCIDWettingPhase){
// 			boundary = true;
// 			  }
// 
// 			if (cell->vertex(ngb)->info().id() == BCIDNonWettingPhase){
// 			boundary = true;
// 			}
 		}
       
       
       if (boundary == false){
       Vt = Vt + max(0.0,(1.0 / cell->info().invVoidVolume())); //std::abs((( std::abs(cell->info().volume()) - std::abs(solver->volumeSolidPore(cell) ) )));
       if (cell->info().saturation() > 0.0){
	Vw = Vw + max(0.0,(1.0*cell->info().saturation() / cell->info().invVoidVolume())); //std::abs((( std::abs(cell->info().volume()) - std::abs(solver->volumeSolidPore(cell) ) )))*cell->info().saturation();
       }
     }
     }
 return (Vw/Vt);  
}

void SoluteFlowEngine::copySaturation()
{	//Copy Saturation list to temporary saturation list to avoid twice movement of a advancing front
    FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles)
    {
	cell->info().saturation() = cell->info().saturationTEMP();
    }
}

void SoluteFlowEngine::copySaturationTEMP()
{	//Copy temporary saturation list back to saturation list

    FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles)
    {
      	cell->info().saturationTEMP() = cell->info().saturation();
      	if(cell->info().saturation() != cell->info().saturationTEMP()){activeiterator = true;}

    }
}


void SoluteFlowEngine::getInsideEqSphere()
{
  
    // This routine finds the radius of the inscribed sphere within each pore-body
    // Following Mackay et al., 1972. 
    // NOTE (thomas): Has to become more efficient
      double d01 = 0.0, d02 = 0.0, d03 = 0.0, d12 = 0.0, d13 = 0.0, d23 = 0.0, Rin = 0.0, r0 = 0.0, r1 = 0.0, r2 =0.0, r3 = 0.0;
      bool check = false;
      unsigned int i = 0;
      double summ = 0.0, count = 0.0;
      
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
	Rin  = 0.0;
	check = false;
	//cout << endl << cell -> info().id;
	//Iterate untill with increasing inscribed sphere radius, untill determinant is zero.
	while (check == false){
	Rin = 0.0 + (min(r0,min(r1,min(r2,r3))) / 1000.0)*i;
	i = i + 1;
	
	M(4,0) = pow((r0+Rin),2);
	M(4,1) = pow((r1+Rin),2);
	M(4,2) = pow((r2+Rin),2);
	M(4,3) = pow((r3+Rin),2);
	M(0,4) = pow((r0+Rin),2);
	M(1,4) = pow((r1+Rin),2);
	M(2,4) = pow((r2+Rin),2);
	M(3,4) = pow((r3+Rin),2);
	
	if (M.determinant() < 0.0){check = true;} //Check whether determinant is negative, if yes, then stop iteration
	if (Rin > 100.0*min(r0,min(r1,min(r2,r3)))){ // Check for upper limits NOTE if this is used, for boundary cells
	  check = true;
	  cell -> info().error()=true;
	  Rin = (cell -> neighbor(0) -> info().InsideSphereRadius() + 
		cell -> neighbor(1) -> info().InsideSphereRadius() +
		cell -> neighbor(2) -> info().InsideSphereRadius() +
		cell -> neighbor(3) -> info().InsideSphereRadius()) / 4.0;
	  //cout << endl << "error with pore " << cell -> info().id;
	}
	}
	
	
	cell -> info().InsideSphereRadius() = Rin;


      }
    //Assign for each pore-body with an error, for the average of surrounding pores. Only applies for B.C pores.  
    FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles){
      if(cell->info().error() == true){ 
    summ = 0.0;
    count = 0;
    if(cell -> neighbor(0) -> info().error() == false){summ = summ + cell -> neighbor(0) -> info().InsideSphereRadius(); count = count+1.0;}
    if(cell -> neighbor(1) -> info().error() == false){summ = summ + cell -> neighbor(1) -> info().InsideSphereRadius(); count = count+1.0;}
    if(cell -> neighbor(2) -> info().error() == false){summ = summ + cell -> neighbor(2) -> info().InsideSphereRadius(); count = count+1.0;}
    if(cell -> neighbor(3) -> info().error() == false){summ = summ + cell -> neighbor(3) -> info().InsideSphereRadius(); count = count+1.0;}
    cell -> info().InsideSphereRadius() = summ / count;
    if (cell -> info().InsideSphereRadius() == 0.0){cell -> info().InsideSphereRadius();}
     
     
   }
  }

      
      
  
    }
    
  

// void SoluteFlowEngine::LiquidBridge()
// {
//     FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles){
//       for (unsigned i = 0; i < 6; i++){
// 	if(cell->info().liquidbridge() [i]){ cerr << "tralalal";}
// 	
//       }
// }
// }
// 

// double SoluteFlowEngine:: liquidbridgeCalc(double R1, double R2, double Pc){
// 	
// 	vector<double> data(2);
// 	double R3 = 0.0;
// 	double tetha1 = 0.0;
// 	double tetha2 = 0.0;
// 	double tetha3 = 0.0;
// 	double surfacetension = 7.2;
// 	double pi = 3.14159265359;
// 	double A = 0.0, B = 0.0, C = 0.0, D = 0.0, E = 0.0;
// 
// 	// Calculate tetha1 from R3
// 	R3 = 2 * surfacetension / (Pc);
// 	tetha1 = acos(R3*(R2 - R1) - R1*(R1 + R2) / (-R3*(R1 + R2) - R1*(R1 + R2)));
// 	if (tetha1 > (0.5*pi)){ tetha1 = 2*pi - tetha1; }
// 	//Calculate all other angles
// 	tetha2 = 2 * atan((R1 / R2)*tan(tetha1 / 2.0));
// 	tetha3 = pi - tetha1 - tetha2;
// 
// 
// 	//Calculate A-E as in Rose (1958)
// 		A = pow((R2 + R3), 2)*(R1 + R2)*(pow(sin(tetha2), 2)) / 6.0;
// 		B = (pow((2 * R3), 2)*(sin(tetha3 / 2.0))*sin((tetha3 / 2.0) + tetha2)) / 3.0;
// 		C = (2.0 / 3.0)*pow(R1, 3)*pow(sin(tetha1 / 2.0), 2);
// 		D = (2.0 / 3.0)*pow(R2, 3)*pow(sin(tetha2 / 2.0), 2);
// 		E = (tetha3 / 2.0)*(R2 + R3)*sin(tetha2)*pow(R3, 2);
// 
// 		data[0] = 2 * pi*(A + B - C - D - E);
// // 		data[1] = 2 * pi*tetha3*R3*((R2 + R3)*sin(tetha2) - (R3*sin((tetha3 / 2.0) + tetha2)*(sin(tetha3 / 2.0) / (tetha3 / 2.0))));
// 
// 	
// 
// 	return data[0];
// 
// }
// 



YADE_PLUGIN ( ( SoluteFlowEngine ) );

#endif //SOLUTE_FLOW
#endif //FLOW_ENGINE

#endif /* YADE_CGAL */
