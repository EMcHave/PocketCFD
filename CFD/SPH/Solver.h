#pragma once
#include "SParticle.h"
#include "HeighbourSearch.h"
namespace SPH
{
	class Solver
	{
	public:
		float2 g = float2(0, -gravity);
		int N = 1000; /*particles per 1 m^2*/
		float stiffness = 1.f;
		float density = 100.f;
		float viscosity = 3.5;
		int kernel_particles = 20;
		float cr = 0.f;
		std::vector<float> Velocities;
		unsigned int counter = 0;

		std::unique_ptr<Geometry> Area;

		Solver(int L, int H, int N, float st, float den, float vis, int kerr_p, float dt, float g);

		void TimeStep();

		std::vector<std::shared_ptr<SParticle>>& Particles() { return m_particles; }

	private:
		
		float particleMass = density / N;
		float particleRadius = std::sqrt(particleMass / M_PI / density);
		float supportRadius = std::sqrt(kernel_particles * particleMass / M_PI / density);
		float smoothing = std::sqrt(kernel_particles * particleMass / M_PI / density);

		float selfDensity = particleMass * (4.f / (M_PI * pow(smoothing, 8))) * pow(smoothing, 6);

		std::vector<std::shared_ptr<SParticle>> m_particles;
		std::vector<std::shared_ptr<SBody>> m_bodies;
		bool m_isFirstStep;

		float m_dt = pow(10, -2);


		void initMidv();
		bool isCollision(std::shared_ptr<SParticle> p);
		float sqDistPoint(std::shared_ptr<SParticle> p);
		void EvaluateForces(uint32_t* table);
		uint32_t* EvaluateHashes();

		float2 ClosestPoint(std::shared_ptr<SParticle> p);

		inline float DefaultKernel(float2 R, float h)
		{
			float r = float2::abs(R);
			if (r <= h)
				return 4.f / (M_PI * pow(h, 8)) * pow(h * h - r * r, 3);
			else
				return 0.f;
		}

		inline float2 GradientKernel(float2 R, float h, float r)
		{		
			if (r <= h)
				return  -30.f / (M_PI * pow(h, 5)) * R / r * pow(h - r, 2);
			else
				return float2(0.f,0.f);
		}

		inline float LaplasianKernel(float2 R, float h, float r)
		{
			if (r <= h)
				return 40.f / (M_PI * pow(h, 5)) * (h - r);
			else 
				return 0.f;
		}
	};
}


