// 2008 © Václav Šmilauer <eudoxos@arcig.cz> 
#pragma once
#include<core/PartialEngine.hpp>
#include<core/State.hpp>
#include<core/Scene.hpp>

class StepDisplacer: public PartialEngine {
	public:
		virtual void action() {
			FOREACH(Body::id_t id, ids){
			const shared_ptr<Body>& b=Body::byId(id,scene);
			if(setVelocities){
				const Real& dt=scene->dt;
				b->state->vel=mov/dt;
				AngleAxisr aa(rot); aa.axis().normalize();
				b->state->angVel=aa.axis()*aa.angle()/dt;
				LOG_DEBUG("Angular velocity set to "<<aa.axis()*aa.angle()/dt<<". Axis="<<aa.axis()<<", angle="<<aa.angle());
			}
			if(!setVelocities){
				b->state->pos+=mov;
				b->state->ori=rot*b->state->ori;
			}
		}
	}
	YADE_CLASS_BASE_DOC_ATTRS(StepDisplacer,PartialEngine,"Apply generalized displacement (displacement or rotation) stepwise on subscribed bodies. Could be used for purposes of contact law tests (by moving one sphere compared to another), but in this case, see rather :yref:`LawTester`",
		((Vector3r,mov,Vector3r::Zero(),,"Linear displacement step to be applied per iteration, by addition to :yref:`State.pos`."))
		((Quaternionr,rot,Quaternionr::Identity(),,"Rotation step to be applied per iteration (via rotation composition with :yref:`State.ori`)."))
		((bool,setVelocities,false,,"If false, positions and orientations are directly updated, without changing the speeds of concerned bodies. If true, only velocity and angularVelocity are modified. In this second case :yref:`integrator<NewtonIntegrator>` is supposed to be used, so that, thanks to this Engine, the bodies will have the prescribed jump over one iteration (dt)."))
	);
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(StepDisplacer);
