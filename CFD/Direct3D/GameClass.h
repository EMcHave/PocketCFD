#pragma once

#include <memory>
#include "Common/DirectXSample.h"
#include "Common/DeviceResources.h"
#include "Common/StepTimer.h"
#include "LogicClass.h"
#include "MoveLookController.h"
#include "Camera.h"
#include "MyRenderer.h"


class GameClass : public winrt::implements<GameClass, winrt::Windows::Foundation::IInspectable>, DX::IDeviceNotify
{
public:
	GameClass(std::shared_ptr<DX::DeviceResources> const& deviceResources);
	~GameClass();
	void CreateWindowSizeDependentResources();
	void CreateDeviceDependentResources();
	void AddObject(std::shared_ptr<Cube>);
	void Update(int i);
	bool Render();

	const std::shared_ptr<LogicClass> Logic() { return m_logic; }

	int delta;

	// IDeviceNotify
	virtual void OnDeviceLost();
	virtual void OnDeviceRestored();

private:

	DX::StepTimer										m_timer;

	std::shared_ptr <LogicClass>						m_logic;

	std::shared_ptr<DX::DeviceResources>				m_deviceResources;

	std::shared_ptr<MyRenderer>							m_sceneRenderer;

	std::shared_ptr<MoveLookController>                 m_controller;

};



