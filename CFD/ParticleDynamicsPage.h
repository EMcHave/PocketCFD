#pragma once

#include "pch.h"
#include "Direct3D/GameClass.h"
#include "ParticleDynamicsPage.g.h"

#include <chrono>
#include <thread>

using namespace winrt::Windows::System::Threading;
using namespace std::chrono_literals;


namespace winrt::CFD::implementation
{
    using namespace Microsoft::Graphics::Canvas::Geometry;
    struct ParticleDynamicsPage : ParticleDynamicsPageT<ParticleDynamicsPage>
    {
    private:
        bool                                        m_isAnimating;
        bool                                        m_timeStepEnded;
        //winrt::agile_ref<CoreWindow> m_window;
        std::shared_ptr<DX::DeviceResources>        m_deviceResources;
        winrt::com_ptr<GameClass>                   m_main;

        winrt::com_ptr<ID3D11VertexShader>          vertexShader;
        winrt::com_ptr<ID3D11PixelShader>           pixelShader;
        winrt::com_ptr<ID3D11InputLayout>           inputLayout;

        winrt::com_ptr<ID3D11Buffer>                vertexBuffer;
        winrt::com_ptr<ID3D11Buffer>                indexBuffer;

        winrt::Windows::Foundation::IAsyncAction    _renderLoopWorker;

        std::vector<float> U_plot, K_plot, E_plot;
        float minY = 0;
        float maxY = 0;
        float dY = 0;
        float resScale = 0;
        float plotWidth;
        float plotHeight;

        CanvasGeometry Path_U{ nullptr };
        CanvasGeometry Path_K{ nullptr };
        CanvasGeometry Path_E{ nullptr };

        
    public:
        ParticleDynamicsPage(){}
        void swapChainPanel_Loaded(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::RoutedEventArgs const& e);
        void PlotCanvas_Draw(winrt::Microsoft::Graphics::Canvas::UI::Xaml::ICanvasAnimatedControl const& sender, winrt::Microsoft::Graphics::Canvas::UI::Xaml::CanvasAnimatedDrawEventArgs const& args);
        void PlotCanvas_Update(winrt::Microsoft::Graphics::Canvas::UI::Xaml::ICanvasAnimatedControl const& sender, winrt::Microsoft::Graphics::Canvas::UI::Xaml::CanvasAnimatedUpdateEventArgs const& args);
        void PlotCanvas_Loaded(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::RoutedEventArgs const& e);
        void PlotCanvas_SizeChanged(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::SizeChangedEventArgs const& e);
    };
}

namespace winrt::CFD::factory_implementation
{
    struct ParticleDynamicsPage : ParticleDynamicsPageT<ParticleDynamicsPage, implementation::ParticleDynamicsPage>
    {
    };
}
