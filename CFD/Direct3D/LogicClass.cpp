#include "pch.h"
#include "LogicClass.h"


#define pidiv4 Directx::XM_PIDIV4;

using namespace DirectX;
LogicClass::LogicClass()
{

}

void LogicClass::Init(
	std::shared_ptr<MoveLookController> const& controller,
	std::shared_ptr<MyRenderer> const& renderer)
{
	m_controller = controller;
	m_renderer = renderer;
	m_camera.SetProjParams(70.0f * 3.1415 / 180.0f, 1.0f, 0.01f, 100.0f);
	m_gravity = XMFLOAT3(0, -0.1, 0);
	m_identityVector = XMFLOAT3(1, 1, 1);

	E = 2 * pow(10, 4);
	G = 2 * pow(10, 2);
	nu = 0.3;
	float rho = 900;
	m = 1;
	a = 0.5;
	J = 1;
	Jp = 5 * pow(10, -10);
	S = 0.0025;
	k = pow(10, 8);



	int NX = 10;
	int NY = 20;
	int NZ = 1;

	InitialConditions initCond1 { XMFLOAT3(0, 0, 8), XMFLOAT3(0.2, -0, 0) };
	InitialConditions initCond2 { XMFLOAT3(10, 0, 0), XMFLOAT3(0, 0, 0) };
	InitialConditions initCond3 { XMFLOAT3(6, 0, 6), XMFLOAT3(0, 0, -0.5) };
	Material mat{ E, G, nu, rho };
	
	auto entity1 = std::make_shared<Entity>(XMFLOAT4(NX, NY, NZ, a), initCond1, mat, EntityType::PlaneTriangle);
	//auto entity2 = std::make_shared<Entity>(XMFLOAT4(NX, NY, NZ, a), initCond3, mat, EntityType::PlaneTriangle);
	auto entity3 = std::make_shared<Entity>(XMFLOAT4(1.5, 11, 10, a), initCond2, mat, EntityType::Cylinder);


	m_entities.push_back(entity1);
	//m_entities.push_back(entity2);
	m_entities.push_back(entity3);
	for (auto entity : m_entities)
	{
		for (int i = 0; i < entity->Particles().size(); i++)
			m_particles.push_back(entity->Particles()[i]);
		
		for (int i = 0; i < entity->Connections().size(); i++)
			m_connections.push_back(entity->Connections()[i]);
			
	}



	m_isFirstTimeStep = true;
}

void LogicClass::TimeStep()
{
	using namespace DirectX;
	
	for (auto entity : m_entities)
	{
		for (auto it = entity->Connections().begin(); it!=entity->Connections().end(); )
		{
			auto connection = *it;
			float dist = GeometricObject::DistanceIJ(connection->p1(), connection->p2());
			if (abs((a - dist)) / a > entity->criticalDeformation)
			{
				connection->Break(entity->ContactParticles());
				it = entity->Connections().erase(it);
				continue;
			}

			XMVECTOR r_ij = connection->p2()->VectorPosition() - connection->p1()->VectorPosition();
			XMVECTOR e_ij = r_ij / dist;
			XMVECTOR f_ij = entity->B1() * (dist - a) * e_ij +
				entity->B2() / 2 / dist * (connection->n_i2(1) - connection->n_i1(1) - XMVector3Dot(e_ij, connection->n_i2(1) - connection->n_i1(1)) * e_ij);
			XMVECTOR MTB = entity->B3() * XMVector3Cross(connection->n_i2(1), connection->n_i1(1)) -
				entity->B4() * 0.5 * (XMVector3Cross(connection->n_i2(2), connection->n_i1(2)) + XMVector3Cross(connection->n_i2(3), connection->n_i1(3)));
			XMVECTOR m_ij = -entity->B2() * 0.5 * XMVector3Cross(e_ij, connection->n_i1(1)) + MTB;
			XMVECTOR m_ji = entity->B2() * 0.5 * XMVector3Cross(e_ij, connection->n_i2(1)) - MTB;

			connection->p1()->Force() += f_ij;
			connection->p2()->Force() -= f_ij;
			connection->p1()->Moment() += m_ij;
			connection->p2()->Moment() += m_ji;

			++it;
			
			//connection->p1()->Moment() += XMVector3Transform(m_ij, XMMatrixTranspose(XMMatrixRotationQuaternion(connection->p1()->VectorQuaternion())));
			//connection->p2()->Moment() += XMVector3Transform(m_ji, XMMatrixTranspose(XMMatrixRotationQuaternion(connection->p2()->VectorQuaternion())));
			
		}

		for (int i = 0; i < (int)entity->Size().x; i++)
			for (int j = 0; j < (int)entity->Size().y; j++)
				if (entity->Particle(i, j, 0)->Position().y < 0)
				{
					auto arch_force = XMFLOAT3(0,
						-5 * m_gravity.y * abs(entity->Particle(i, j, 0)->Position().y * (entity->Size().x - 1) * (entity->Size().y - 1) * entity->Size().w),
						0);
					entity->Particle(i, j, 0)->Force() += XMLoadFloat3(&arch_force);
				}


		for (auto another_entity : m_entities)
			if (another_entity != entity)
				for (auto particle : entity->ContactParticles())
					for (auto another_particle : another_entity->ContactParticles())
					{
						float r = GeometricObject::DistanceIJ(particle, another_particle);
						float f = Entity::LennardJones(r, a, 0.005);
						DirectX::XMVECTOR e_ij = (another_particle->VectorPosition() - particle->VectorPosition()) / r;
						particle->Force() += f * e_ij;
					}

		
		for (auto particle : entity->Particles())
		{
			particle->Force() += XMLoadFloat3(&m_gravity);

			if (particle->Position().y > 0)
				vis = 0;
			else
				vis = 0;

			XMFLOAT3 constraint(1, 1, 1);
			if (m_isFirstTimeStep)
			{
				particle->Velocity(particle->VectorVelocity() + 0.5 * dt * particle->Force() / m);
				//particle->AngularVelocity(particle->VectorVelocity() + 0.5 * dt * particle->Moment() / J);
			}
				

			TranslationalIntegratorLeapFrog(particle);
			RotationalIntegratorSymplectic(particle);
			
			if (!m_realTimeRender)
				particle->AddTimePoint(particle->VectorPosition(), particle->VectorQuaternion());
		}
		
		for (auto connection : entity->Connections())
			connection->Update();
		for (auto particle : entity->Particles())
			particle->ResetForcesAndMoments();
		if (m_isFirstTimeStep)
			m_isFirstTimeStep = false;
	}
	
}

void LogicClass::AnimationStep(int i)
{
	for (auto p : m_particles)
	{
		p->Position(p->m_positions[i]);
		p->Quaternion(p->m_orientations[i]);
	}
}

void LogicClass::EvaluateEnergy(
	std::vector<float>& U,
	std::vector<float>& K,
	std::vector<float>& E)
{
	float U_total = 0.0f;
	float K_total = 0.0f;

	DirectX::XMFLOAT3 tempf3;
	float tempf;

	auto e = m_entities[0];

	for (auto p : e->Particles())
	{
		auto v2om2 = XMVector3Dot(p->VectorVelocity(), p->VectorVelocity()) + 
			XMVector3Dot(p->VectorAngularVelocity(), p->VectorAngularVelocity());
		XMStoreFloat3(&tempf3, v2om2);
		K_total += 0.5 * tempf3.x;
	}
	for (auto it = e->Connections().begin(); it != e->Connections().end(); ++it)
	{
		auto c = *it;
		if (c.get() != nullptr)
		{
			float dist = GeometricObject::DistanceIJ(c->p2(), c->p1());
			XMVECTOR r_ij = c->p2()->VectorPosition() - c->p1()->VectorPosition();
			XMVECTOR e_ij = r_ij / dist;

			XMVECTOR dU = 0.5 * e->B1() * (dist - e->Size().w) * (dist - e->Size().w) * DirectX::XMLoadFloat3(&m_identityVector) +
				0.5 * e->B2() * XMVector3Dot(c->v2_i(0)->VectorOrientation() - c->v1_i(0)->VectorOrientation(), e_ij) +
				e->B3() * XMVector3Dot(c->v2_i(0)->VectorOrientation(), c->v1_i(0)->VectorOrientation()) -
				0.5 * e->B4() * (DirectX::XMVector3Dot(c->v2_i(1)->VectorOrientation(), c->v1_i(1)->VectorOrientation()) + DirectX::XMVector3Dot(c->v2_i(2)->VectorOrientation(), c->v1_i(2)->VectorOrientation()));
			XMStoreFloat(&tempf, dU);

			U_total += tempf + e->B2() + e->B3() + e->B4();
		}
	}
	
	U.push_back(U_total);
	K.push_back(K_total);
	E.push_back(U_total + K_total);
}

void LogicClass::RotationalIntegratorSymplectic(const std::shared_ptr<Cube> particle)
{
	
	XMVECTOR av_prev = particle->VectorAngularVelocity() * !m_isFirstTimeStep + particle->Moment() / J * (0.5 * dt);
	XMStoreFloat3(&m_tempAngularVelocity3, av_prev);
	m_tempAngularVelocity4 = XMFLOAT4(m_tempAngularVelocity3.x, m_tempAngularVelocity3.y, m_tempAngularVelocity3.z, 0.f);
	XMVECTOR q_half_plus = particle->VectorQuaternion() + (0.5 * dt) * 0.5 * XMQuaternionMultiply(particle->VectorQuaternion(), XMLoadFloat4(&m_tempAngularVelocity4));

	XMVECTOR av_half_plus = particle->VectorAngularVelocity() * !m_isFirstTimeStep + dt * particle->Moment() / J;
	XMStoreFloat3(&m_tempAngularVelocity3, av_half_plus);
	m_tempAngularVelocity4 = XMFLOAT4(m_tempAngularVelocity3.x, m_tempAngularVelocity3.y, m_tempAngularVelocity3.z, 0.f);
	XMVECTOR q_new = particle->VectorQuaternion() + 0.5 * dt * XMQuaternionMultiply(q_half_plus, XMLoadFloat4(&m_tempAngularVelocity4));

	particle->AngularVelocity(av_half_plus);
	particle->Quaternion(q_new / XMQuaternionLength(q_new));
	
	/*
	XMVECTOR L_prev_global = particle->KineticMoment() + (0.5 * dt) * particle->Moment();
	XMVECTOR L_prev_local = XMVector3Transform(L_prev_global, XMMatrixRotationQuaternion(particle->VectorQuaternion()));
	XMStoreFloat3(&m_tempAngularVelocity3, L_prev_local/J);
	m_tempAngularVelocity4 = XMFLOAT4(m_tempAngularVelocity3.x, m_tempAngularVelocity3.y, m_tempAngularVelocity3.z, 0.f);
	XMVECTOR q_half_plus = particle->VectorQuaternion() + (0.5 * dt) * 0.5 * XMQuaternionMultiply(particle->VectorQuaternion(), XMLoadFloat4(&m_tempAngularVelocity4));

	XMVECTOR L_half_plus_global = particle->KineticMoment() + dt * particle->Moment();
	XMVECTOR L_half_plus_local = XMVector3Transform(L_half_plus_global, XMMatrixRotationQuaternion(q_half_plus));
	XMStoreFloat3(&m_tempAngularVelocity3, L_half_plus_local / J);
	m_tempAngularVelocity4 = XMFLOAT4(m_tempAngularVelocity3.x, m_tempAngularVelocity3.y, m_tempAngularVelocity3.z, 0.f);
	XMVECTOR q_new = particle->VectorQuaternion() + 0.5 * dt * XMQuaternionMultiply(q_half_plus, XMLoadFloat4(&m_tempAngularVelocity4));

	particle->KineticMoment(L_half_plus_global);
	particle->Quaternion(q_new / XMQuaternionLength(q_new));
	*/
}

void LogicClass::RotationalIntegratorNonSymplectic(const std::shared_ptr<Cube> particle)
{
	XMVECTOR av_new = particle->VectorAngularVelocity() + particle->Moment() / J * dt;
	XMStoreFloat3(&m_tempAngularVelocity3, av_new);
	m_tempAngularVelocity4 = XMFLOAT4(m_tempAngularVelocity3.x, m_tempAngularVelocity3.y, m_tempAngularVelocity3.z, 0.f);
	XMVECTOR q_new = particle->VectorQuaternion() + 0.5 * XMQuaternionMultiply(XMLoadFloat4(&m_tempAngularVelocity4), particle->VectorQuaternion()) * dt;

	particle->AngularVelocity(av_new);
	particle->Quaternion(q_new / XMQuaternionLength(q_new));
}

void LogicClass::TranslationalIntegratorLeapFrog(const std::shared_ptr<Cube> particle)
{
	XMVECTOR v_new = particle->VectorVelocity() + (-vis * particle->VectorVelocity() + (!particle->m_isFixed) * particle->Force() / m) * dt * !m_isFirstTimeStep;
	XMVECTOR r_new = particle->VectorPosition() + v_new * dt;

	particle->Velocity(v_new);
	particle->Position(r_new);
}

