#include "pch.h"
#include "Solver.h"

namespace SPH
{
	Solver::Solver(int L, int H, int N, float st, float den, float vis, int kerr_p, float dt, float g)
		: N(N), stiffness(st), density(den), viscosity(vis), kernel_particles(kerr_p), m_dt(dt)
	{
		this->g = SPH::float2(0, -g);

		Area = std::make_unique<Geometry>(L, H);
		m_bodies.push_back(
			std::make_shared<SBody>(particleMass, density, 1.f, 2.f, float2(0.f,0.f), float2(0.f, 0.f))
		);

		Velocities.reserve(m_particles.size());

		for (auto body : m_bodies)
		{
			for (auto p : body->m_particles)
				m_particles.push_back(p);
		}
		initMidv();
	}

	void Solver::initMidv()
	{
		uint32_t* table = EvaluateHashes();
		EvaluateForces(table);
		for (auto& p : m_particles)
			p->mid_v(-0.5 * m_dt * p->f() / p->rho());
	}

	void Solver::TimeStep()
	{
		uint32_t* table = EvaluateHashes();
		Velocities.clear();

		EvaluateForces(table);
		for (auto& p : m_particles)
		{
			float2 acceleration = p->f() / p->rho();
			p->v(p->v() + acceleration * m_dt);
			p->r(p->r() + p->v() * m_dt);
			Velocities.push_back(float2::abs(p->v()));
			/*
			float2 new_mid_v = p->mid_v() + m_dt * (acceleration);
			float2 new_r = p->r() + m_dt * new_mid_v;

			p->r(new_r);
			p->mid_v(new_mid_v);
			p->v((new_mid_v + p->mid_v()) / 2);
			*/


			
			if (isCollision(p))
			{
				float2 cp(0.f,0.f);
				float2 n(0.f,0.f);
				cp = ClosestPoint(p);
				//float pen = float2::abs(p->r() - cp);

				if (cp.x == 0 || cp.x == Area->L) { n.x = (cp.x == 0) ? -1 : 1; }
				if (cp.y == 0 || cp.y == Area->H) { n.y = (cp.y == 0) ? -1 : 1; }

				p->r(cp);
				p->v(p->v() - p->v().dot(n) * n);
			}	
			
		}
		counter++;	
	}

	void Solver::EvaluateForces(uint32_t* table)
	{

		for (size_t pIndex = 0; pIndex < m_particles.size(); pIndex++)
		{
			float rho = 0.f;
			std::shared_ptr<SParticle> p = m_particles[pIndex];
			p->ResetForce();
			int2 cell = getCell(p, smoothing);

			for (int x = -1; x <= 1; x++)
			{
				for (int y = -1; y <= 1; y++)
				{
					uint32_t cellHash = getHash(cell + int2(x, y));
					uint32_t apIndex = table[cellHash];

					if (apIndex == NO_PARTICLE) {
						continue;
					}
					while (apIndex < m_particles.size()) {
						if (pIndex == apIndex) {
							apIndex++;
							continue;
						}
						std::shared_ptr<SParticle> ap = m_particles[apIndex];
						if (ap->hash != cellHash) {
							break;
						}
						rho += particleMass * DefaultKernel(p->r() - ap->r(), smoothing);
						apIndex++;
					}
				}
			}

			p->rho(selfDensity + rho);
			p->p(stiffness * (p->rho() - density));
		}

		for (size_t pIndex = 0; pIndex < m_particles.size(); pIndex++)
		{
			std::shared_ptr<SParticle> p = m_particles[pIndex];
			float2 f_pressure;
			float2 f_viscosity;
			int2 cell = getCell(p, smoothing);

			for (int x = -1; x <= 1; x++)
			{
				for (int y = -1; y <= 1; y++)
				{
					uint32_t cellHash = getHash(cell + int2(x, y));
					uint32_t apIndex = table[cellHash];

					if (apIndex == NO_PARTICLE) {
						continue;
					}
					while (apIndex < m_particles.size()) {
						if (pIndex == apIndex) {
							apIndex++;
							continue;
						}
						std::shared_ptr<SParticle> ap = m_particles[apIndex];
						if (ap->hash != cellHash) {
							break;
						}
						float dist = (float2::abs(p->r() - ap->r()));
						if(dist > pow(10, -7))
						{
							f_pressure += -(p->p() + ap->p()) / 2 * (particleMass / ap->rho()) * GradientKernel(p->r() - ap->r(), smoothing, dist);
							f_viscosity += viscosity * (ap->v() - p->v()) * (particleMass / ap->rho()) * LaplasianKernel(ap->r() - p->r(), smoothing, dist);
						}	
						apIndex++;
					}
				}
			}

			float2 totalForce = g * p->rho() + f_pressure + f_viscosity;

			p->f(totalForce);
		}

		delete(table);
	}

	bool Solver::isCollision(std::shared_ptr<SParticle> p)
	{
		int fx[2]{ 0,0 };
		for (int i = 0; i < 2; i++)
		{
			float v = p->r()[i];
			fx[i] = (v <= Area->min[i]) || (v >= Area->max[i]);
		}
		return fx[0] || fx[1];
	}

	float Solver::sqDistPoint(std::shared_ptr<SParticle> p)
	{
		float sqDist = 0.0f;

		for (int i = 0; i < 2; i++)
		{
			float v = (i == 0) ? p->r().x : p->r().y;

			if (v < Area->min[i]) sqDist += (Area->min[i] - v) * (Area->min[i] - v);
			if (v > Area->max[i]) sqDist += (v - Area->max[i]) * (v - Area->max[i]);
		}

		return sqDist;
	}

	

	uint32_t* Solver::EvaluateHashes()
	{
		for (auto& p : m_particles)
			p->hash = getHash(getCell(p, smoothing));

		std::sort(std::execution::par_unseq,
			begin(m_particles), end(m_particles),
			[&](const std::shared_ptr<SParticle>& i, const std::shared_ptr<SParticle>& j) {
				return i->hash < j->hash;
			}
		);

		return createNeighborTable(m_particles);
	}

	float2 Solver::ClosestPoint(std::shared_ptr<SParticle> p)
	{
		float2 point(0, 0);

		for (int i = 0; i < 2; i++)
		{
			float v = p->r()[i];

			if (v <= Area->min[i]) v = Area->min[i];
			if (v >= Area->max[i]) v = Area->max[i];
			point[i] = v;
		}
		return point;
	}
}

