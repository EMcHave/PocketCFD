#pragma once

#include "SPHPage.g.h"
#include "SPH/Solver.h"

using namespace winrt::Windows::Foundation::Numerics;

namespace winrt::CFD::implementation
{
    struct SPHPage : SPHPageT<SPHPage>
    {
        SPHPage()
        {

        }
        void Solver() {}
    public:
        void solveButton_Click(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::RoutedEventArgs const& e);
    private:
        std::unique_ptr<SPH::Solver> m_solver;
        double x_scale;
        double y_scale;
        double height;
        double width;
        float maxV = 0;

        const float3 RainBow7[7] = {MY_INDIGO, MY_BLUE, MY_LIGHT_BLUE, MY_GREEN, MY_YELLOW, MY_ORANGE, MY_RED };
    public:
        winrt::Windows::UI::Color ValueColorBessonov(double x, int IsStriped);
        void canvas_Draw(winrt::Microsoft::Graphics::Canvas::UI::Xaml::ICanvasAnimatedControl const& sender, winrt::Microsoft::Graphics::Canvas::UI::Xaml::CanvasAnimatedDrawEventArgs const& args);
        void canvas_Update(winrt::Microsoft::Graphics::Canvas::UI::Xaml::ICanvasAnimatedControl const& sender, winrt::Microsoft::Graphics::Canvas::UI::Xaml::CanvasAnimatedUpdateEventArgs const& args);
    };
}

namespace winrt::CFD::factory_implementation
{
    struct SPHPage : SPHPageT<SPHPage, implementation::SPHPage>
    {
    };
}
