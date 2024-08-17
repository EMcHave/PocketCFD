#include "pch.h"
#include "GameClass.h"


GameClass::GameClass(std::shared_ptr<DX::DeviceResources> const& deviceResources)
	:m_deviceResources(deviceResources)
{
	m_deviceResources->RegisterDeviceNotify(this);
	m_sceneRenderer = std::make_shared<MyRenderer>(m_deviceResources);
	m_controller = std::make_shared<MoveLookController>(winrt::Windows::UI::Core::CoreWindow::GetForCurrentThread());
	m_logic = std::make_shared<LogicClass>();
	m_logic->Init(m_controller, m_sceneRenderer);
	delta = 0;
	m_timer.SetFixedTimeStep(true);
	m_timer.SetTargetElapsedSeconds(1.f / 60.f);
}

GameClass::~GameClass()
{
	m_deviceResources->RegisterDeviceNotify(nullptr);
}

void GameClass::CreateWindowSizeDependentResources()
{
	m_sceneRenderer->CreateWindowSizeDependentResources();
}

void GameClass::CreateDeviceDependentResources()
{
	m_sceneRenderer->CreateDeviceDependentResources(m_logic);
}

void GameClass::AddObject(std::shared_ptr<Cube> cube)
{
	//m_objects.push_back(cube);
}

void GameClass::Update(int i)
{
	m_timer.Tick([&]()
		{
			m_controller->Update();
			m_logic->GameCamera().SetViewParams(
				m_controller->Position(),
				m_controller->Orientation(),
				DirectX::XMFLOAT3(0, 1, 0)
			);
			if (m_logic->IsRealTime())
				m_logic->TimeStep();
			else
				m_logic->AnimationStep(i);
		});
}

bool GameClass::Render()
{
	m_sceneRenderer->Render();

	return true;
}

void GameClass::OnDeviceLost()
{
	m_sceneRenderer->ReleaseDeviceDependentResources();
}

void GameClass::OnDeviceRestored()
{
	m_sceneRenderer->CreateDeviceDependentResources(m_logic);
	CreateWindowSizeDependentResources();
}
