#include "pch.h"
#include "SPHPage.h"
#if __has_include("SPHPage.g.cpp")
#include "SPHPage.g.cpp"
#endif

using namespace winrt;
using namespace Windows::UI::Xaml;
using namespace Windows::UI;
using namespace winrt::Microsoft::UI::Xaml::Controls;
using namespace Microsoft::Graphics::Canvas::UI::Xaml;
using namespace Windows::Foundation;
using namespace Windows::Foundation::Numerics;

namespace winrt::CFD::implementation
{

}


void winrt::CFD::implementation::SPHPage::solveButton_Click(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::RoutedEventArgs const& e)
{
    maxV = 0;
    m_solver = std::make_unique<SPH::Solver>(
        lField().Value(), hField().Value(), static_cast<int>(NField().Value()),
        stField().Value(), rhoField().Value(), nuField().Value(),
        static_cast<int>(kerpField().Value()), dtField().Value(), gField().Value());
    canvas().Paused(true);
    float ratio = m_solver->Area->L / m_solver->Area->H;
    if (ratio > 1)
    {
        canvas().Width(canvasBorder().ActualWidth() * 0.9);
        canvas().Height(canvas().Width() * 0.9 / ratio);
        height = canvas().Height();
        width = canvas().Width();
    }
    else
    {
        canvas().Height(canvasBorder().ActualHeight() * 0.9);
        canvas().Width(canvasBorder().ActualWidth() * 0.9 * ratio);
        height = canvas().Height();
        width = canvas().Width();
    }


    canvas().Visibility(Visibility::Visible);

    x_scale = width / m_solver->Area->L;
    y_scale = height / m_solver->Area->H;

    canvas().Paused(false);
}


void winrt::CFD::implementation::SPHPage::canvas_Draw(winrt::Microsoft::Graphics::Canvas::UI::Xaml::ICanvasAnimatedControl const& sender, winrt::Microsoft::Graphics::Canvas::UI::Xaml::CanvasAnimatedDrawEventArgs const& args)
{
    args.DrawingSession().DrawLine(5, 5, width - 5,  5, Colors::Green());
    args.DrawingSession().DrawLine(5, 5, 5, height - 5, Colors::Green());
    args.DrawingSession().DrawLine(width - 5, 5, width - 5, height - 5, Colors::Green());
    args.DrawingSession().DrawLine(5, height -  5, width -  5, height - 5, Colors::Green());
    
    for (auto p : m_solver->Particles())
    {
     
        Color col = ValueColorBessonov(SPH::float2::abs(p->v()) / maxV, 0);
        args.DrawingSession().FillCircle(p->r().x * x_scale, height - (p->r().y) * y_scale, p->R * x_scale, col);
    }
}



void winrt::CFD::implementation::SPHPage::canvas_Update(winrt::Microsoft::Graphics::Canvas::UI::Xaml::ICanvasAnimatedControl const& sender, winrt::Microsoft::Graphics::Canvas::UI::Xaml::CanvasAnimatedUpdateEventArgs const& args)
{
    m_solver->TimeStep();
    std::sort(std::execution::par_unseq, begin(m_solver->Velocities), end(m_solver->Velocities));
    float new_maxV = m_solver->Velocities.back();
    if (new_maxV > maxV) { maxV = new_maxV; }
}

winrt::Windows::UI::Color winrt::CFD::implementation::SPHPage::ValueColorBessonov(double x, int IsStriped)
{

    if (x <= 0) x = 0.001;
    else if (x >= 1) x = 0.999;
    double N7 = 7; 
    double Y = x * (N7 - 1);
    int J = (int)Y;
    float3 RRR_;
    if (IsStriped)
        RRR_ = RainBow7[J] * 255;
    else
        RRR_ = (RainBow7[J] + (RainBow7[J + 1] - RainBow7[J]) * (Y - J)) * 255;
    return ColorHelper::FromArgb(255, (uint8_t)RRR_.x, (uint8_t)RRR_.y, (uint8_t)RRR_.z);
}
