#include "pch.h"
#include "ParticleDynamicsPage.h"
#if __has_include("ParticleDynamicsPage.g.cpp")
#include "ParticleDynamicsPage.g.cpp"
#endif

namespace winrt::CFD::implementation
{


}

void winrt::CFD::implementation::ParticleDynamicsPage::swapChainPanel_Loaded(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::RoutedEventArgs const& e)
{
    {
        using namespace winrt::Windows::System::Threading;
        using namespace winrt::Windows::Foundation;

        auto scale = swapChainPanel().CompositionScaleX();

        m_deviceResources = std::make_shared<DX::DeviceResources>();
        m_deviceResources->SetSwapChain(swapChainPanel());
        if (!m_main)
            m_main = make_self<GameClass>(m_deviceResources);

        m_main->CreateDeviceDependentResources();
        m_main->CreateWindowSizeDependentResources();

        m_main->Logic()->IsRealTime(1);
        m_main->Logic()->DT(pow(10, -1));

        m_isAnimating = true;
        m_timeStepEnded = false;

        auto workItemHandler = WorkItemHandler([this](IAsyncAction action)
            {
                if (m_main->Logic()->IsRealTime())
                {
                    while (true)
                    {
                        if (m_isAnimating)
                        {
                            m_timeStepEnded = false;
                            m_main->Update(0);
                            if (m_main->Render())
                                m_deviceResources->Present();
                            m_timeStepEnded = true;
                        }
                    }
                }
                else
                {
                    float time = 2;
                    int N = time / m_main->Logic()->DT();

                    int dt = N / (time * 60);

                    for (int n = 0; n < N; n++)
                        m_main->Logic()->TimeStep();

                    while (true)
                    {
                        for (int n = 0; n < N; n += dt)
                        {
                            m_main->Update(n);
                            if (m_main->Render())
                                m_deviceResources->Present();
                        }
                    }
                }
            }
        );

        _renderLoopWorker = ThreadPool::RunAsync(workItemHandler, WorkItemPriority::High, WorkItemOptions::TimeSliced);
    }
}

void winrt::CFD::implementation::ParticleDynamicsPage::PlotCanvas_Draw(winrt::Microsoft::Graphics::Canvas::UI::Xaml::ICanvasAnimatedControl const& sender, winrt::Microsoft::Graphics::Canvas::UI::Xaml::CanvasAnimatedDrawEventArgs const& args)
{
    args.DrawingSession().DrawGeometry(Path_U, winrt::Windows::UI::Colors::Red());
    args.DrawingSession().DrawGeometry(Path_K, winrt::Windows::UI::Colors::Green());
    args.DrawingSession().DrawGeometry(Path_E, winrt::Windows::UI::Colors::Cyan());
    args.DrawingSession().DrawText(L"U - energy", float2(plotWidth - 150, 10), winrt::Windows::UI::Colors::Red());
    args.DrawingSession().DrawText(L"K - energy", float2(plotWidth - 150, 40), winrt::Windows::UI::Colors::Green());
    args.DrawingSession().DrawText(L"Total energy", float2(plotWidth - 150, 70), winrt::Windows::UI::Colors::Cyan());
    //args.DrawingSession().DrawCircle(float2(0, 0), 40, winrt::Windows::UI::Colors::Red());
}


void winrt::CFD::implementation::ParticleDynamicsPage::PlotCanvas_Update(winrt::Microsoft::Graphics::Canvas::UI::Xaml::ICanvasAnimatedControl const& sender, winrt::Microsoft::Graphics::Canvas::UI::Xaml::CanvasAnimatedUpdateEventArgs const& args)
{
    using namespace winrt::Microsoft::Graphics::Canvas;
    m_isAnimating = false;
    while (!m_timeStepEnded){}
        //std::this_thread::sleep_for(0.1ms);
    m_main->Logic()->EvaluateEnergy(U_plot, K_plot, E_plot);
    float max = max(U_plot.back(), max(K_plot.back(), E_plot.back()));
    if (max > maxY) 
        maxY = max;

    float min = min(U_plot.back(), min(K_plot.back(), E_plot.back()));
    if (min < minY) 
        minY = min;

    dY = abs(maxY - minY);
    resScale = plotHeight / dY;

    Geometry::CanvasPathBuilder build_U(PlotCanvas().as<ICanvasResourceCreator>());
    Geometry::CanvasPathBuilder build_K(PlotCanvas().as<ICanvasResourceCreator>());
    Geometry::CanvasPathBuilder build_E(PlotCanvas().as<ICanvasResourceCreator>());

    build_U.BeginFigure(0, plotHeight - resScale * (U_plot[0] - minY));
    build_K.BeginFigure(0, plotHeight - resScale * (K_plot[0] - minY));
    build_E.BeginFigure(0, plotHeight - resScale * (E_plot[0] - minY));

    for (int it = 1; it < U_plot.size(); it++)
    {
        double u = plotHeight - resScale * (U_plot[it] - minY);
        double k = plotHeight - resScale * (K_plot[it] - minY);
        double e = plotHeight - resScale * (E_plot[it] - minY);
        build_U.AddLine(it * plotWidth / U_plot.size(), u);
        build_K.AddLine(it * plotWidth / K_plot.size(), k);
        build_E.AddLine(it * plotWidth / E_plot.size(), e);
    }

    build_U.EndFigure(Geometry::CanvasFigureLoop::Open);
    build_K.EndFigure(Geometry::CanvasFigureLoop::Open);
    build_E.EndFigure(Geometry::CanvasFigureLoop::Open);
    Path_U = Geometry::CanvasGeometry::CreatePath(build_U);
    Path_K = Geometry::CanvasGeometry::CreatePath(build_K);
    Path_E = Geometry::CanvasGeometry::CreatePath(build_E);
    m_isAnimating = true;
}


void winrt::CFD::implementation::ParticleDynamicsPage::PlotCanvas_Loaded(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::RoutedEventArgs const& e)
{
    PlotCanvas().Paused(false);
}


void winrt::CFD::implementation::ParticleDynamicsPage::PlotCanvas_SizeChanged(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::SizeChangedEventArgs const& e)
{
    plotHeight = PlotCanvas().ActualHeight();
    plotWidth = PlotCanvas().ActualWidth();
}
