#pragma once
#include <ppl.h>
#include "Camera.h"
#include "Cube.h"
#include "Line.h"
#include "Entity.h"
#include "MoveLookController.h"
#include "Connection.h"


class MyRenderer;



class LogicClass
{
public:
	LogicClass();

	void Init(
		_In_ std::shared_ptr<MoveLookController> const& controller,
		_In_ std::shared_ptr<MyRenderer> const& renderer
	);
	Camera& GameCamera() { return m_camera; };
	std::vector<std::shared_ptr<Cube>> const& RenderObjects() { return m_particles; };
	std::vector<std::shared_ptr<Connection>> const& RenderVectors() { return m_connections; };
	std::shared_ptr<Cube>& Object(int i) { return m_particles.at(i); };

	void TimeStep();
	void AnimationStep(int i);
	void EvaluateEnergy(
		std::vector<float>& U,
		std::vector<float>& K,
		std::vector<float>& E);

	void DT(float dt) { this->dt = dt; }
	float DT() { return dt; }

	bool IsRealTime() { return m_realTimeRender; }
	void IsRealTime(bool isRealTime) { m_realTimeRender = isRealTime; }

private:
	
	void RotationalIntegratorSymplectic(const std::shared_ptr<Cube> particle);
	void RotationalIntegratorNonSymplectic(const std::shared_ptr<Cube> particle);
	void TranslationalIntegratorLeapFrog(const std::shared_ptr<Cube> particle);

	std::shared_ptr<MoveLookController>				m_controller;
	std::shared_ptr<MyRenderer>						m_renderer;
	Camera											m_camera;
	std::vector<std::shared_ptr<Cube>>				m_particles;
	std::vector<std::shared_ptr<Entity>>			m_entities;
	std::vector<std::shared_ptr<Connection>>		m_connections;

	XMFLOAT3										m_gravity;
	XMFLOAT3										m_identityVector;
	XMFLOAT4										m_tempAngularVelocity4;
	XMFLOAT3										m_tempAngularVelocity3;
	
	bool											m_realTimeRender;
	bool											m_isFirstTimeStep;

	float dt;
	float vis;

	float m;
	float E;
	float nu;
	float J;
	float Jp;
	float G;
	float S;
	float a;
	float k;
};

