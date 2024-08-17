#include "pch.h"
#include "Connection.h"

Connection::Connection(
	std::shared_ptr<Cube> p1,
	std::shared_ptr<Cube> p2,
	std::array<std::shared_ptr<Line>, 3> v1,
	std::array<std::shared_ptr<Line>, 3> v2)
	: m_particle1(p1),
	m_particle2(p2),
	m_vectors1(v1),
	m_vectors2(v2)
{
	m_isBroken = false;
	p1->RegisterConnection();
	p2->RegisterConnection();
}

void Connection::Update()
{
	if (!m_isBroken)
	{
		for (std::shared_ptr<Line> vec : m_vectors1)
		{
			vec->Quaternion(p1()->VectorQuaternion());
			vec->Position(p1()->VectorPosition());
		}
		for (std::shared_ptr<Line> vec : m_vectors2)
		{
			vec->Quaternion(p2()->VectorQuaternion());
			vec->Position(p2()->VectorPosition());
		}
	}
}

void Connection::Break(std::vector<std::shared_ptr<Cube>>& contactParticles)
{
	m_isBroken = true;
	m_particle1->UnRegisterConnection();
	m_particle2->UnRegisterConnection();
	if (!m_particle1->IsContact())
	{
		contactParticles.push_back(m_particle1);
		m_particle1->MakeContact();
	}
	if (!m_particle2->IsContact())
	{
		contactParticles.push_back(m_particle2);
		m_particle2->MakeContact();
	}
}
