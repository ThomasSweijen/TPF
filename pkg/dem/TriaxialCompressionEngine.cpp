/*************************************************************************
*  Copyright (C) 2006 by Bruno Chareyre                                  *
*  bruno.chareyre@hmg.inpg.fr                                            *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "TriaxialCompressionEngine.hpp"
#include<core/Scene.hpp>
#include<core/Omega.hpp>
#include<lib/base/Math.hpp>
#include<boost/lambda/lambda.hpp>
#include<pkg/dem/Shop.hpp>
#include<core/Interaction.hpp>
#include<pkg/common/Sphere.hpp>
#include<pkg/dem/FrictPhys.hpp>
#include<pkg/common/ElastMat.hpp>

class Ip2_CohFrictMat_CohFrictMat_CohFrictPhys;

CREATE_LOGGER(TriaxialCompressionEngine);
YADE_PLUGIN((TriaxialCompressionEngine));

TriaxialCompressionEngine::~TriaxialCompressionEngine()
{
}

void TriaxialCompressionEngine::doStateTransition(stateNum nextState){
	if ( /* currentState==STATE_UNINITIALIZED && */ nextState==STATE_ISO_COMPACTION){
		sigma_iso=sigmaIsoCompaction;
		previousSigmaIso=sigma_iso;
	}
	else if(nextState==STATE_TRIAX_LOADING){
		sigma_iso=sigmaLateralConfinement;
		previousSigmaIso=sigma_iso;
		internalCompaction = false;
		if (frictionAngleDegree>0) setContactProperties(frictionAngleDegree);
		height0 = height; depth0 = depth; width0 = width;
		//compressionActivated = true;
		wall_bottom_activated=false; wall_top_activated=false;
		if(currentState==STATE_ISO_UNLOADING && !noFiles){ LOG_INFO("Speres -> /tmp/unloaded.spheres"); Shop::saveSpheresToFile("/tmp/unloaded.spheres"); }
		if(!firstRun && !noFiles) saveSimulation=true; // saving snapshot .xml will actually be done in ::action
		Phase1End = "Unloaded";
	}
	else if(currentState==STATE_ISO_COMPACTION && nextState==STATE_ISO_UNLOADING){
		sigma_iso=sigmaLateralConfinement;
		sigmaIsoCompaction = sigmaLateralConfinement;
		previousSigmaIso=sigma_iso;
		internalCompaction=false; // unloading will not change grain sizes
		if (frictionAngleDegree>0) setContactProperties(frictionAngleDegree);
		if(!firstRun && !noFiles) saveSimulation=true;
		Phase1End = "Compacted";
	}
	else if ((currentState==STATE_ISO_COMPACTION || currentState==STATE_ISO_UNLOADING) && nextState==STATE_LIMBO) {
	//urrentState==STATE_DIE_COMPACTION
		internalCompaction = false;
		if (frictionAngleDegree>0) setContactProperties(frictionAngleDegree);
		height0 = height; depth0 = depth; width0 = width;
		if(!noFiles) saveSimulation=true; // saving snapshot .xml will actually be done in ::action
		// stop simulation here, since nothing will happen from now on
		Phase1End = (currentState==STATE_ISO_COMPACTION ? "compacted" : "unloaded");
		if(!noFiles) Shop::saveSpheresToFile("/tmp/limbo.spheres");
		// Please keep this saving process intact, I'm tired of running 3 days simulations and getting nothing at the end!
		if(!firstRun && !noFiles) saveSimulation=true; // saving snapshot .xml will actually be done in ::action
	}
	else if( nextState==STATE_FIXED_POROSITY_COMPACTION){
		internalCompaction = false;
		wall_bottom_activated=false; wall_top_activated=false;
		wall_front_activated=false; wall_back_activated=false;
		wall_right_activated=false; wall_left_activated=false;
	}
	else { LOG_ERROR("Undefined transition from "<<stateName(currentState)<<" to "<<stateName(nextState)<<"! (ignored)"); return; }

	LOG_INFO("State transition from "<<stateName(currentState)<<" to "<<stateName(nextState)<<" done.");
	currentState=nextState;
	previousState=currentState; // should be always kept in sync, used to track manual changes to the .xml
}

void TriaxialCompressionEngine::updateParameters ()
{
	UnbalancedForce=ComputeUnbalancedForce ();

	if (  (currentState!=STATE_TRIAX_LOADING && currentState==STATE_ISO_COMPACTION) || currentState==STATE_ISO_UNLOADING || currentState==STATE_FIXED_POROSITY_COMPACTION || autoCompressionActivation)
	{
		if (UnbalancedForce<=StabilityCriterion && std::abs ( ( meanStress-sigma_iso ) /sigma_iso ) <0.005 && fixedPoroCompaction==false )
		{
			// only go to UNLOADING if it is needed
			if ( currentState==STATE_ISO_COMPACTION && autoUnload && sigmaLateralConfinement!=sigmaIsoCompaction ) {
				doStateTransition (STATE_ISO_UNLOADING );
				computeStressStrain (); // update stress and strain
			}
			else if((currentState==STATE_ISO_COMPACTION || currentState==STATE_ISO_UNLOADING || currentState==STATE_LIMBO) && autoCompressionActivation){
				doStateTransition (STATE_TRIAX_LOADING );
				computeStressStrain (); // update stress and strain
			}
		}
		else if ( porosity<=fixedPorosity && currentState==STATE_FIXED_POROSITY_COMPACTION )
		{
// 			Omega::instance().pause();
			return;
		}
	}
}

void TriaxialCompressionEngine::action()
{
	if (!warn++) LOG_WARN ("This engine is deprecated, please switch to TriaxialStressController if you expect long term support.")
	// here, we make sure to get consistent parameters, in case someone fiddled with the scene .xml manually
	if ( firstRun )
	{
		LOG_INFO ( "First run, will initialize!" );
		if ( (sigmaIsoCompaction!=previousSigmaIso || currentState==STATE_UNINITIALIZED || currentState== STATE_LIMBO) && currentState!=STATE_TRIAX_LOADING && !fixedPoroCompaction) doStateTransition (STATE_ISO_COMPACTION );
		if (previousState!=STATE_TRIAX_LOADING && currentState==STATE_TRIAX_LOADING)
			doStateTransition (STATE_TRIAX_LOADING );
		if (fixedPorosity<1 && currentState==STATE_UNINITIALIZED && fixedPoroCompaction) doStateTransition (STATE_FIXED_POROSITY_COMPACTION );
		previousState=currentState;
		previousSigmaIso=sigma_iso;
		firstRun=false; // change this only _after_ state transitions
	}
	if ( scene->iter % testEquilibriumInterval == 0 )
	{
		updateParameters ();
		maxStress = max(maxStress,stress[wall_top][1]);
		LOG_INFO("UnbalancedForce="<< UnbalancedForce<<", rel stress "<< std::abs ( ( meanStress-sigma_iso ) /sigma_iso ));
	}
	if ( saveSimulation )
	{
		if(!noFiles){
			string fileName = "./"+ Key + "_" + Phase1End + "_" +
							  boost::lexical_cast<string> ( scene->iter ) + "_" +
							  boost::lexical_cast<string> ( currentState ) + ".xml";
			LOG_INFO ( "saving snapshot: "<<fileName );
			Omega::instance().saveSimulation ( fileName );
			fileName="./"+ Key + "_"+Phase1End+"_"+boost::lexical_cast<string> ( scene->iter ) + "_" +
					 boost::lexical_cast<string> ( currentState ) +".spheres";
			LOG_INFO ( "saving spheres: "<<fileName );
			Shop::saveSpheresToFile ( fileName );
		}
		saveSimulation = false;
	}
	if (isAxisymetric || internalCompaction){
		if (stressMask & 1) goal1=sigma_iso;
		if (stressMask & 2) goal2=sigma_iso;
		if (stressMask & 3) goal3=sigma_iso;
	}
	
	TriaxialStressController::action();
	if ( currentState==STATE_LIMBO && autoStopSimulation )
	{
		Omega::instance().pause();
		return;
	}


	if ( currentState==STATE_TRIAX_LOADING )
	{
		if ( scene->iter % 100 == 0 )
		{
			LOG_INFO ("Triax Compression started");
		}
		if (scene->iter % 100 == 0) LOG_DEBUG("Compression active.");
		const Real& dt = scene->dt;

		if (std::abs(epsilonMax) > std::abs(strain[1])) {
			if ( currentStrainRate != strainRate ) currentStrainRate += ( strainRate-currentStrainRate ) *0.0003;
			/* Move top and bottom wall according to strain rate */
			State* p_bottom=Body::byId(wall_bottom_id,scene)->state.get();
			p_bottom->pos += 0.5*currentStrainRate*height*translationAxis*dt;
			State* p_top=Body::byId(wall_top_id,scene)->state.get();
			p_top->pos -= 0.5*currentStrainRate*height*translationAxis*dt;
		}
	}
	if ( currentState==STATE_FIXED_POROSITY_COMPACTION )
	{
		if ( scene->iter % 100 == 0 ) LOG_INFO ("Compression started");
		const Real& dt = scene->dt;
		State* p_bottom=Body::byId(wall_bottom_id,scene)->state.get();
		State* p_top=Body::byId(wall_top_id,scene)->state.get();
		State* p_left=Body::byId(wall_left_id,scene)->state.get();
		State* p_right=Body::byId(wall_right_id,scene)->state.get();
		State* p_front=Body::byId(wall_front_id,scene)->state.get();
		State* p_back=Body::byId(wall_back_id,scene)->state.get();

		/* Move top and bottom wall according to strain rate */
		p_bottom->pos += 0.5*strainRate*height*translationAxis*dt;
		p_top->pos -= 0.5*strainRate*height*translationAxis*dt;
		p_back->pos += 0.5*strainRate*depth*translationAxisz*dt;
		p_front->pos -= 0.5*strainRate*depth*translationAxisz*dt;
		p_left->pos += 0.5*strainRate*width*translationAxisx*dt;
		p_right->pos -= 0.5*strainRate*width*translationAxisx*dt;
	}

}

void TriaxialCompressionEngine::setContactProperties(Real frictionDegree)
{
	Shop::setContactFriction(frictionDegree*Mathr::PI/180.0);
}

