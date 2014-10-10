// 2004 © Janek Kozicki <cosurgi@berlios.de> 
// 2009 © Václav Šmilauer <eudoxos@arcig.cz> 
// 2014 © Raphael Maurin <raphael.maurin@irstea.fr> 

#include"ForceEngine.hpp"
#include<yade/core/Scene.hpp>
#include<yade/pkg/common/Sphere.hpp>
#include<yade/lib/smoothing/LinearInterpolate.hpp>
#include<yade/pkg/dem/Shop.hpp>

#include<yade/core/IGeom.hpp>
#include<yade/core/IPhys.hpp>

#include <boost/random/linear_congruential.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

YADE_PLUGIN((ForceEngine)(InterpolatingDirectedForceEngine)(RadialForceEngine)(DragEngine)(LinearDragEngine)(HydroForceEngine));

void ForceEngine::action(){
	FOREACH(Body::id_t id, ids){
		if (!(scene->bodies->exists(id))) continue;
		scene->forces.addForce(id,force);
	}
}

void InterpolatingDirectedForceEngine::action(){
	Real virtTime=wrap ? Shop::periodicWrap(scene->time,*times.begin(),*times.rbegin()) : scene->time;
	direction.normalize(); 
	force=linearInterpolate<Real,Real>(virtTime,times,magnitudes,_pos)*direction;
	ForceEngine::action();
}

void RadialForceEngine::postLoad(RadialForceEngine&){ axisDir.normalize(); }

void RadialForceEngine::action(){
	FOREACH(Body::id_t id, ids){
		if (!(scene->bodies->exists(id))) continue;
		const Vector3r& pos=Body::byId(id,scene)->state->pos;
		Vector3r radial=(pos - (axisPt+axisDir * /* t */ ((pos-axisPt).dot(axisDir)))).normalized();
		if(radial.squaredNorm()==0) continue;
		scene->forces.addForce(id,fNorm*radial);
	}
}

void DragEngine::action(){
	FOREACH(Body::id_t id, ids){
		Body* b=Body::byId(id,scene).get();
		if (!b) continue;
		if (!(scene->bodies->exists(id))) continue;
		const Sphere* sphere = dynamic_cast<Sphere*>(b->shape.get());
		if (sphere){
			Real A = sphere->radius*sphere->radius*Mathr::PI;	//Crossection of the sphere
			Vector3r velSphTemp = b->state->vel;
			Vector3r dragForce = Vector3r::Zero();
			
			if (velSphTemp != Vector3r::Zero()) {
				dragForce = -0.5*Rho*A*Cd*velSphTemp.squaredNorm()*velSphTemp.normalized();
			}
			scene->forces.addForce(id,dragForce);
		}
	}
}

void LinearDragEngine::action(){
	FOREACH(Body::id_t id, ids){
		Body* b=Body::byId(id,scene).get();
		if (!b) continue;
		if (!(scene->bodies->exists(id))) continue;
		const Sphere* sphere = dynamic_cast<Sphere*>(b->shape.get());
		if (sphere){
			Vector3r velSphTemp = b->state->vel;
			Vector3r dragForce = Vector3r::Zero();
			
			Real b = 6*Mathr::PI*nu*sphere->radius;
			
			if (velSphTemp != Vector3r::Zero()) {
				dragForce = -b*velSphTemp;
			}
			scene->forces.addForce(id,dragForce);
		}
	}
}

void HydroForceEngine::action(){
	/* Velocity fluctuation determination (not usually done at each dt, that is why it not placed in the other loop) */
	if (velFluct == true){
		/* check size */
		size_t size=vFluct.size();
		if(size<scene->bodies->size()){
			size=scene->bodies->size();
			vFluct.resize(size);
		}
		/* reset stored values to zero */
		memset(& vFluct[0],0,sizeof(Vector2r)*size);
	
		/* Create a random number generator rnd() with a gaussian distribution of mean 0 and stdev 1.0 */
		/* see http://www.boost.org/doc/libs/1_55_0/doc/html/boost_random/reference.html and the chapter 7 of Numerical Recipes in C, second edition (1992) for more details */
		static boost::minstd_rand0 randGen((int)TimingInfo::getNow(true));
		static boost::normal_distribution<Real> dist(0.0, 1.0);
		static boost::variate_generator<boost::minstd_rand0&,boost::normal_distribution<Real> > rnd(randGen,dist);

		double rand1 = 0.0;
		double rand2 = 0.0;
		/* Attribute a fluid velocity fluctuation to each body above the bed elevation */
		FOREACH(Body::id_t id, ids){
			Body* b=Body::byId(id,scene).get();
			if (!b) continue;
			if (!(scene->bodies->exists(id))) continue;
			const Sphere* sphere = dynamic_cast<Sphere*>(b->shape.get());
			if (sphere){
				Vector3r posSphere = b->state->pos;//position vector of the sphere
				int p = floor((posSphere[2]-zRef)/deltaZ); //cell number in which the particle is
				// if the particle is inside the water and above the bed elevation, so inside the turbulent flow, evaluate a turbulent fluid velocity fluctuation which will be used to apply the drag.
				if ((p<nCell)&&(posSphere[2]-zRef>bedElevation)) { 
					Real uStar2 = simplifiedReynoldStresses[p];
					if (uStar2>0.0){
						Real uStar = sqrt(uStar2);
						rand1 = rnd();
						rand2 = -rand1 + rnd();
						vFluct[id] = Vector2r(rand1*uStar,rand2*uStar);
					}
				}
			}
			
		}
		velFluct = false;
	}
	
	/* Application of hydrodynamical forces */
	FOREACH(Body::id_t id, ids){
		Body* b=Body::byId(id,scene).get();
		if (!b) continue;
		if (!(scene->bodies->exists(id))) continue;
		const Sphere* sphere = dynamic_cast<Sphere*>(b->shape.get());
		if (sphere){
			Vector3r posSphere = b->state->pos;//position vector of the sphere
			int p = floor((posSphere[2]-zRef)/deltaZ); //cell number in which the particle is
			if ((p<nCell)&&(p>0)) {
				Vector3r liftForce = Vector3r::Zero();
				Vector3r dragForce = Vector3r::Zero();
				Vector2r fluctVelBody = vFluct[id];//fluid velocity fluctuation associated to the particle's position considered.
				Vector3r vFluid(vxFluid[p]+fluctVelBody.x(),0.0,fluctVelBody.y()); //fluid velocity at the particle's position
				Vector3r vPart = b->state->vel;//particle velocity
				Vector3r vRel = vFluid - vPart;//fluid-particle relative velocity

				//Drag force calculation
				Real Rep = vRel.norm()*sphere->radius*2*rhoFluid/viscoDyn; //particles reynolds number
				Real A = sphere->radius*sphere->radius*Mathr::PI;	//Crossection of the sphere
				if (vRel.norm()!=0.0) {
					Real hindranceF = pow(1-phiPart[p],-expoRZ); //hindrance function
					Real Cd = (0.44 + 24.4/Rep)*hindranceF; //drag coefficient
					dragForce = 0.5*rhoFluid*A*Cd*vRel.squaredNorm()*vRel.normalized();
				}
				//lift force calculation due to difference of fluid pressure between top and bottom of the particle
				int intRadius = floor(sphere->radius/deltaZ);
				if ((p+intRadius<nCell)&&(p-intRadius>0)&&(lift==true)) {
					Real vRelTop = vxFluid[p+intRadius] - vPart[0]; // relative velocity of the fluid wrt the particle at the top of the particle
					Real vRelBottom = vxFluid[p-intRadius] - vPart[0]; // same at the bottom
					liftForce[2] = 0.5*rhoFluid*A*Cl*(vRelTop*vRelTop-vRelBottom*vRelBottom);
				}
				//buoyant weight force calculation
				Vector3r buoyantForce = -4.0/3.0*Mathr::PI*sphere->radius*sphere->radius*sphere->radius*rhoFluid*gravity;
				//add the hydro forces to the particle
				scene->forces.addForce(id,dragForce+liftForce+buoyantForce);		
			}
		}
	}
}

