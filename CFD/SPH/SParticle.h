#pragma once
#include "includes.h"

namespace SPH
{
	struct float2
	{
		float x;
		float y;

		float2() : x(0), y(0){}
		float2(float x, float y);

		float min_val() { return min(x, y); }
		float max_val() { return max(x, y); }

		static float abs(float2 const& r1);
		float dot(float2& f2);



		float2 operator+ (float2) const;
		float2 operator- (float2) const;

		float2& operator+= (float2);
		float2& operator-= (float2);

		float& operator[](int i);

		friend float2 operator* (float2, float);
		friend float2 operator* (float, float2);

		friend float2 operator/ (float2, float);
		friend float2 operator/ (float, float2);
	};

	struct int2
	{
		int x;
		int y;

		int2() : x(0), y(0) {}
		int2(int x, int y) : x(x), y(y) {}

		int2 operator+ (int2 i2) const { return int2(x + i2.x, y + i2.y); }
		int2 operator- (int2 i2) const { return int2(x - i2.x, y - i2.y); }
	};

	struct Geometry
	{
		float L;
		float H;

		const float min[2]{ 0.f, 0.f };
		const float max[2]{ L, H };

		Geometry() : L(1.f), H(1.f) {}
		Geometry(float l, float h): L(l), H(h) {}

		float2 Center() { return float2(L / 2, H / 2); }
	};

	struct SParticle
	{
		float R;

		SParticle() {}
		SParticle(float2 r, float2 v, float R);
		uint32_t hash;
		float2 r() { return m_position; }
		float2 v() { return m_velocity; }
		float2 mid_v() { return m_midvelocity; }
		//float2 u() { return m_displacement; }
		float rho() { return m_density; }
		float p() { return m_pressure; }

		float2& f() { return m_force; }

		void r(float2 pos) { m_position = pos; }
		void v(float2 vel) { m_velocity = vel; }
		void mid_v(float2 vel) { m_midvelocity = vel; }
		//void u(float2 disp) { m_displacement = disp; }
		void rho(float rho) { m_density = rho; }
		void p(float p) { m_pressure = p; }
		void f(float2 f) { m_force = f; }

		void ResetForce() { m_force = float2(0, 0); }

	private:
		unsigned int id;

		float2 m_position;
		float2 m_midvelocity;
		float2 m_velocity;
		//float2 m_displacement;

		float m_density;
		float m_pressure;

		float2 m_force;
	};

	struct SBody
	{
		unsigned int id;
		const float s = 0.8;
		float m_particle_radius;
		std::vector<std::shared_ptr<SParticle>> m_particles;

		SBody(
			float m,
			float rho,
			float lx,
			float ly,
			float2 pos,
			float2 vel
		) 
		{
			m_particle_radius = std::sqrt(m/(M_PI * rho));
			CreateBlock(int(lx / (2 * m_particle_radius)), int(ly /(2 * m_particle_radius)), pos, vel);
		}
	
		void CreateBlock(unsigned int nx, unsigned int ny, float2 pos, float2 vel)
		{


			for (int j  = 0; j < ny; j++)
				for (int i = 0; i < nx; i++)
				{
					m_particles.push_back(
						std::make_shared<SParticle>(
							float2(pos.x + m_particle_radius + i * 2 * m_particle_radius, pos.y + m_particle_radius + j * 2 * m_particle_radius),
							float2(vel.x, vel.y),
							s * m_particle_radius
						)
					);
				}
		}
	};

}


