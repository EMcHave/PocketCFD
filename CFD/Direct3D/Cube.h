#pragma once

#include "Line.h"
using namespace DirectX;
class Cube : public GeometricObject
{
public:
	Cube();
	Cube(DirectX::XMFLOAT3 pos, float, bool);
	Cube(DirectX::XMFLOAT3 pos, DirectX::XMVECTOR quat, float, bool);

	void SetInitialQuaternion(XMVECTOR);

	using GeometricObject::Position;
	void Position(DirectX::XMFLOAT3 position);
	void Position(DirectX::XMVECTOR position);
	void Quaternion(DirectX::XMFLOAT4 quaternion);
	void Quaternion(DirectX::XMVECTOR quaternion);
	DirectX::XMFLOAT4 Quaternion();

	void RegisterConnection() { m_numberOfConnections += 1; }
	void UnRegisterConnection() { m_numberOfConnections -= 1; }
	uint32_t NumberOfConnections() { return m_numberOfConnections; }

	void Force(DirectX::XMVECTOR f) { m_force = f; }
	DirectX::XMVECTOR& Force() { return m_force; }

	void Moment(DirectX::XMVECTOR m) { m_moment = m; }
	DirectX::XMVECTOR& Moment() { return m_moment; }

	bool IsContact() { return m_isContact; }
	void MakeContact() { m_isContact = true; }

	void ResetForcesAndMoments()
	{
		m_force = XMVectorZero();
		m_moment = XMVectorZero();
	}
	
	void AddTimePoint(XMVECTOR pos, XMVECTOR quat) 
	{ 
		m_positions.push_back(pos); 
		m_orientations.push_back(quat);
	}

	std::vector<XMVECTOR>					m_positions;
	std::vector<XMVECTOR>					m_orientations;
private:
	void Update();
	


	DirectX::XMVECTOR						m_force;
	DirectX::XMVECTOR						m_moment;
	uint32_t								m_numberOfConnections;
	float									m_scale;
	bool									m_isContact;
};

