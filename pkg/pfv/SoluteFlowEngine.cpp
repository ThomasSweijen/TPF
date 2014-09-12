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
	Real bubblePressure2;
	bool imbibition2;
	bool networkConnectivity2;
	std::vector<bool> liquidbr;
	
	
	
// 	liquidbr.resize(4, 0);inline std::vector<double>& Rh (void) {return rayHydr;}
 	SoluteCellInfo (void) : FlowCellInfo_SoluteFlowEngineT() {WSaturation2=0.0;InsideSphereRadius2=0.0;bubblePressure2=0.0;liquidbr.resize(6, 0);}
 	inline Real& saturationTEMP (void) {return WSaturation2;}
 	inline const Real& saturationTEMP (void) const {return WSaturation2;}
 	inline std::vector<bool>& liquidbridge (void) {return liquidbr;}
 	
  	inline Real& InsideSphereRadius (void) {return InsideSphereRadius2;}
  	inline const Real& InsideSphereRadius (void) const {return InsideSphereRadius2;}
  	inline bool& imbibition (void) {return imbibition2;} //If imbibition is true, an interface is present within pore body
  	inline const bool& imbibition (void) const {return imbibition2;}
  	inline bool& networkConnectivity (void) {return networkConnectivity2;} //If imbibition is true, an interface is present within pore body
  	inline const bool& networkConnectivity (void) const {return networkConnectivity2;}
  	
 	inline void getInfo (const SoluteCellInfo& otherCellInfo) {FlowCellInfo_SoluteFlowEngineT::getInfo(otherCellInfo); saturationTEMP()=otherCellInfo.saturationTEMP();}
};

typedef TemplateFlowEngine_SoluteFlowEngineT<SoluteCellInfo,FlowVertexInfo_SoluteFlowEngineT> SoluteFlowEngineT;
REGISTER_SERIALIZABLE(SoluteFlowEngineT);
YADE_PLUGIN((SoluteFlowEngineT));

class SoluteFlowEngine : public SoluteFlowEngineT
{
	public :
		void getInsideEqSphere();
		void TPFaction();
		void LiquidBridge();
		double liquidbridgeCalc(double R1, double R2, double Pc);
		void copySaturation();
		void copySaturationTEMP();
		void staticIterator();
		void InitializeScenario();
		void InsertBoundaryConditions();
		void checkConnectivityNetwork();
		int printConnectivityNetwork();
		void initializeImbibition();
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

		,,,
		.def("TPFaction",&SoluteFlowEngine::TPFaction,"Activate Two Phase Flow Engine")
		.def("printConnectivityNetwork",&SoluteFlowEngine::printConnectivityNetwork,"show non-connected pores")
		.def("initializeImbibition",&SoluteFlowEngine::initializeImbibition,"Initialize all imbibition bools to false")
		.def("LiquidBridge",&SoluteFlowEngine::LiquidBridge,"Calculate swelling of liquid bridge")
		.def("getBulkSaturation",&SoluteFlowEngine::getBulkSaturation,"Activate Two Phase Flow Engine")
		.def("staticIterator",&SoluteFlowEngine::staticIterator,"Iterator to find eq. in two-phase flow")
		.def("getSaturation",&SoluteFlowEngine::getSaturation,(boost::python::arg("ID")),"Get Saturation")
		.def("getInscribedRadius",&SoluteFlowEngine::getInscribedRadius,(boost::python::arg("ID")),"Get Saturation")
		.def("getInsideEqSphere",&SoluteFlowEngine::getInsideEqSphere,"Get equivalent pore volume")
		.def("InitializeScenario",&SoluteFlowEngine::InitializeScenario,"Reset Initial conditions")
		.def("InsertBoundaryConditions",&SoluteFlowEngine::InsertBoundaryConditions,"Insert the boundary conditions into vectors.")
		
		)
};
REGISTER_SERIALIZABLE(SoluteFlowEngine);


// PeriodicFlowEngine::~PeriodicFlowEngine(){}
void SoluteFlowEngine::LiquidBridge()
{
    FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles){
      for (unsigned i = 0; i < 6; i++){
	if(cell->info().liquidbridge() [i]){ cerr << "tralalal";}
	
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
    cout << endl << "Calculate radius of inscribed sphere is done!";
    }
  
  active = true;
  
  // While loop for calculating equilibrium conditions for a set Pn-Pw, only solved for continuous flow. 
   while(active){
  checkConnectivityNetwork();
  copySaturation();
  active = false;
  staticIterator();
  InsertBoundaryConditions();
  copySaturationTEMP();
  if(makeMovie){solver->saveVtk("./VTK");}  
  }

  
  
}

void SoluteFlowEngine::InitializeScenario()
{
  FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles)
 	{
 	 cell->info().saturationTEMP() = initialWettingSaturation;
	 cell->info().imbibition() = false;
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
void SoluteFlowEngine::initializeImbibition()
{
    FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles)
    {
	cell->info().imbibition() = false;
    }
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
      	if(cell->info().saturation() != cell->info().saturationTEMP()){active = true;}
	cell->info().saturationTEMP() = cell->info().saturation();
    }
}


void SoluteFlowEngine::staticIterator()
{
	double Rcp = 0.0, throatRadius = 0.0;
	double CapillaryPressure = 0.0;

	CapillaryPressure = BCPressureNonWettingPhase - BCPressureWettingPhase;
	CappilaryPressureOverall = CapillaryPressure;
	Rcp = (2.0*(72.0 / 1000.0)) / CapillaryPressure;
	

	  
	 
	FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles)
	{	  
	  
	  
	  // (1) Check for pore-body imbibition
	  if((cell -> info().imbibition() == true) && (cell->info().networkConnectivity() == true)){
	   if (Rcp >= cell -> info().InsideSphereRadius()){
	    cell -> info().saturation() = 1.0; //NOTE (thomas): might need to include residual NW-phase
	   }
	  }
	  
	  
	  // (2) Check for pore-throat criteria
	  
	  if (cell->info().saturation() != 0.0){ //We only look at water saturated pores
	  
	  for (unsigned int ngb=0;ngb<4;ngb++)
	  {
	  
	    // check for neigboring pore with NW-saturation, and no imbibition interface
	    if ((cell -> neighbor(ngb) ->info().saturation() < 1.0) && (cell->info().imbibition() == false)){
	      if ((cell-> neighbor(ngb) -> info().networkConnectivity() == true) && (cell-> info().networkConnectivity() == true))
	      {
	      
		throatRadius = std::abs(solver->computeHydraulicRadius(cell, ngb));
		//(a) wetting connectivity, or imbibition of neighboring throat
		if (throatRadius <= Rcp){
		  cell -> neighbor(ngb) -> info().saturation() = 0.0; //Although an interface is present in neigboring pore, saturation is zero, residual?
		  cell -> neighbor(ngb) -> info().imbibition() = true;
		}
		//(b) non-wetting connectivity, or drainage of this pore body
		if (throatRadius > Rcp){
		  cell -> info().saturation() = 0.0; //residual water saturation, liquid bridge?
		}
	  
	    }
	    }
	  
	  }
	  }
	  }
} 

int SoluteFlowEngine::printConnectivityNetwork()
{
  int summ =0;
  FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles)
  {
  if(cell->info().networkConnectivity() == true){
    cout << endl << "ID "<< cell->info().id << " Saturation "<< cell->info().saturation();
    summ =summ+1;
  }  
  }
  
  return summ;
}

void SoluteFlowEngine::checkConnectivityNetwork()
{
      bool active = true;
      FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles)
      {
	cell -> info().networkConnectivity() = false;
	for (unsigned int ngb=0;ngb<4;ngb++){ 
	  if (cell->vertex(ngb)->info().id() == BCIDWettingPhase){
	    cell->info().networkConnectivity() = true;
	  }
	  if (cell->vertex(ngb)->info().id() == BCIDNonWettingPhase){
	    cell->info().networkConnectivity() = true;

	  }
	}
      }
      
     
      
      while(active == true){
	active = false;
	
	
      FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles)
      {
	if(cell->info().networkConnectivity() == true){
	  for(unsigned int ngb = 0; ngb < 4; ngb++)
	  {
	    if(cell->neighbor(ngb)->info().saturation() == cell->info().saturation())
	      if ( cell->neighbor(ngb)->info().networkConnectivity() == false){active = true;}
	      cell->neighbor(ngb)->info().networkConnectivity() = true;
	  }
	}
      }

      }
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

       for (unsigned int ngb=0;ngb<4;ngb++){ 
 			if (cell->vertex(ngb)->info().id() == BCIDWettingPhase){
			boundary = true;
			  }

			if (cell->vertex(ngb)->info().id() == BCIDNonWettingPhase){
			boundary = true;
			}
 		}
       
       
       if (boundary == false){
       Vt = Vt + (( std::abs(cell->info().volume()) - std::abs(solver->volumeSolidPore(cell) ) ));
       if (cell->info().saturation() > 0.0){
	Vw = Vw +  ((( std::abs(cell->info().volume()) - std::abs(solver->volumeSolidPore(cell) ) ))*cell->info().saturation());
       }
     }
     }
 return (Vw/Vt);  
}




void SoluteFlowEngine::getInsideEqSphere()
{
  
    // This routine finds the radius of the inscribed sphere within each pore-body
    // Following Mackay et al., 1972. 
    // NOTE (thomas): Has to become more efficient though
      double d01 = 0.0, d02 = 0.0, d03 = 0.0, d12 = 0.0, d13 = 0.0, d23 = 0.0, Rin = 0.0, r0 = 0.0, r1 = 0.0, r2 =0.0, r3 = 0.0;
      bool check = false;
      unsigned int i = 0;
      
      Eigen::MatrixXd M(6,6);


      
      FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles){
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
	

	
	r0 = sqrt(cell -> vertex(0) -> point().weight());
	r1 = sqrt(cell -> vertex(1) -> point().weight());
	r2 = sqrt(cell -> vertex(2) -> point().weight());
	r3 = sqrt(cell -> vertex(3) -> point().weight());
	
	
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
	
	while (check == false){
	Rin = 0.0 + (min(r0,min(r1,min(r2,r3))) / 5000.0)*i;
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
	if (Rin > 10.0*min(r0,min(r1,min(r2,r3)))){ // Check for upper limits NOTE if this is used, an error has occured
	  check = true;
	  Rin = cell -> neighbor(3) -> info().InsideSphereRadius();
	  cout << endl << "error with pore " << cell -> info().id;
	}
	
	
	//cout << endl << "i "<<i  << " pore " << cell -> info().id << " r0 " << r0<< " check " << check << " Rin " << Rin << " Det " << M.determinant();
	}
	
	
	cell -> info().InsideSphereRadius() = Rin;


      }
  }
  
//     double Px = 0.0;
//     double Py = 0.0;
//     double Pz = 0.0;
//     double dx = 0.0;
//     double dy = 0.0;
//     double dz = 0.0;
//     double length1, length2, length3, length4 = 0.0;
//      FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles)
//      {
// 
// 	
//     Px = (cell->vertex(0)->point().x() + cell->vertex(1)->point().x() + cell->vertex(2)->point().x() + cell->vertex(3)->point().x()) / 4.0;
//     Py = (cell->vertex(0)->point().y() + cell->vertex(1)->point().y() + cell->vertex(2)->point().y() + cell->vertex(3)->point().y()) / 4.0;
//     Pz = (cell->vertex(0)->point().z() + cell->vertex(1)->point().z() + cell->vertex(2)->point().z() + cell->vertex(3)->point().z()) / 4.0;
//     
//     for (unsigned int ngb = 0; ngb < 4; ngb ++){
//       dx = cell->vertex(ngb)->point().x() - Px;
//       dy = cell->vertex(ngb)->point().y() - Py;
//       dz = cell->vertex(ngb)->point().z() - Pz;
//       if (cell ->info().id == 10){cerr << endl << dx << endl << dy << endl << dz;}
//       if (ngb == 0){ length1 = sqrt(pow(dx,2) + pow(dy,2) + pow(dz,2)) - sqrt(cell->vertex(ngb)->point().weight());}
//       if (ngb == 1){ length2 = sqrt(pow(dx,2) + pow(dy,2) + pow(dz,2)) - sqrt(cell->vertex(ngb)->point().weight());}
//       if (ngb == 2){ length3 = sqrt(pow(dx,2) + pow(dy,2) + pow(dz,2)) - sqrt(cell->vertex(ngb)->point().weight());}
//       if (ngb == 3){ length4 = sqrt(pow(dx,2) + pow(dy,2) + pow(dz,2)) - sqrt(cell->vertex(ngb)->point().weight());}
//     }
//      
//      
//      cell ->info().InsideSphereRadius() = (length1 + length2 + length3 + length4) / 4.0; 
//      
//          if (cell ->info().id == 10){
// 	   cerr << endl << Px << endl << Py << endl << Pz << endl << length1 << endl << length2 << endl << length3 << endl << length4 << endl << cell ->info().InsideSphereRadius();
// 	 }
//       
//       
//     }
//     
    



double SoluteFlowEngine:: liquidbridgeCalc(double R1, double R2, double Pc){
	
	vector<double> data(2);
	double R3 = 0.0;
	double tetha1 = 0.0;
	double tetha2 = 0.0;
	double tetha3 = 0.0;
	double surfacetension = 7.2;
	double pi = 3.14159265359;
	double A = 0.0, B = 0.0, C = 0.0, D = 0.0, E = 0.0;

	// Calculate tetha1 from R3
	R3 = 2 * surfacetension / (Pc);
	tetha1 = acos(R3*(R2 - R1) - R1*(R1 + R2) / (-R3*(R1 + R2) - R1*(R1 + R2)));
	if (tetha1 > (0.5*pi)){ tetha1 = 2*pi - tetha1; }
	//Calculate all other angles
	tetha2 = 2 * atan((R1 / R2)*tan(tetha1 / 2.0));
	tetha3 = pi - tetha1 - tetha2;


	//Calculate A-E as in Rose (1958)
		A = pow((R2 + R3), 2)*(R1 + R2)*(pow(sin(tetha2), 2)) / 6.0;
		B = (pow((2 * R3), 2)*(sin(tetha3 / 2.0))*sin((tetha3 / 2.0) + tetha2)) / 3.0;
		C = (2.0 / 3.0)*pow(R1, 3)*pow(sin(tetha1 / 2.0), 2);
		D = (2.0 / 3.0)*pow(R2, 3)*pow(sin(tetha2 / 2.0), 2);
		E = (tetha3 / 2.0)*(R2 + R3)*sin(tetha2)*pow(R3, 2);

		data[0] = 2 * pi*(A + B - C - D - E);
// 		data[1] = 2 * pi*tetha3*R3*((R2 + R3)*sin(tetha2) - (R3*sin((tetha3 / 2.0) + tetha2)*(sin(tetha3 / 2.0) / (tetha3 / 2.0))));

	

	return data[0];

}




YADE_PLUGIN ( ( SoluteFlowEngine ) );

#endif //SOLUTE_FLOW
#endif //FLOW_ENGINE

#endif /* YADE_CGAL */
