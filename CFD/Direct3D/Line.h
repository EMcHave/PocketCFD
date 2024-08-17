#pragma once
#include "GeometricObject.h"

class Line : public GeometricObject
{
public:
	Line(Axis axis);
	Line(Axis axis, DirectX::XMFLOAT3 pos, DirectX::XMVECTOR direction);
	Line(Axis axis, DirectX::XMFLOAT3 pos, DirectX::XMVECTOR direction, DirectX::XMVECTOR initQuat);


	DirectX::XMFLOAT3 Orientation() { return m_vector; }
	DirectX::XMVECTOR VectorOrientation() { return DirectX::XMLoadFloat3(&m_vector); }
	void Orientation(DirectX::XMFLOAT3 orient) { m_vector = orient; }
	void Orientation(DirectX::XMVECTOR orient) { XMStoreFloat3(&m_vector, orient); }

	using GeometricObject::Position;
	void Position(DirectX::XMFLOAT3 position);
	void Position(DirectX::XMVECTOR position);

	void Quaternion(DirectX::XMFLOAT4 position);
	void Quaternion(DirectX::XMVECTOR position);

private:
	void Update();
	void Init();
	DirectX::XMFLOAT3 m_vector;
	DirectX::XMFLOAT3 m_initVector;
	Axis m_axis;
};

__forceinline void Line::Position(DirectX::XMFLOAT3 position)
{
	m_position = position;
	Update();
}

__forceinline void Line::Position(DirectX::XMVECTOR position)
{
	DirectX::XMStoreFloat3(&m_position, position);
	Update();
}

__forceinline void Line::Quaternion(DirectX::XMFLOAT4 quaternion)
{
	m_quaternion = quaternion;
	XMStoreFloat3(&m_vector, DirectX::XMVector3Rotate(XMLoadFloat3(&m_vector), XMLoadFloat4(&m_quaternion)));
	Update();
	
}

__forceinline void Line::Quaternion(DirectX::XMVECTOR quaternion)
{
	using namespace DirectX;
	DirectX::XMStoreFloat4(&m_quaternion, quaternion);
	XMStoreFloat3(&m_vector,XMVector3Rotate(XMLoadFloat3(&m_initVector), quaternion));
	Update();
}