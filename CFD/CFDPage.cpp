#include "pch.h"
#include <iostream>
#include "CFDPage.h"
#include <Windows.h>
#include <unknwn.h>
#include <pplawait.h>
#if __has_include("CFDPage.g.cpp")
#include "CFDPage.g.cpp"
#endif

using namespace winrt;
using namespace Windows::UI;
using namespace Windows::UI::Xaml;
using namespace Microsoft::Graphics::Canvas::UI::Xaml;
using namespace winrt::CFD::implementation;

IAsyncAction winrt::CFD::implementation::CFDPage::solveButton_Click(IInspectable const& sender, RoutedEventArgs const& e)
{
    using namespace concurrency;

    progressBar().Visibility(Visibility::Visible());
    debugBuffer += L"Начат расчет\n";
    debugInfo().Text(debugBuffer);
    percentOfSolution().Text(L"0 %");

    int begin = GetTickCount();
    double T = Solver().T();
    solution = Solver().Solve();
    solution.Progress([&](auto const& sender, double progress)
        {
            progressBar().Value(100*progress);
            percentOfSolution().Text(to_wstring((int)(100*progress)) + L" %");
        });
    co_await solution;
    int end = GetTickCount() - begin;
    
    
    percentOfSolution().Text(L"Расчет завершен успешно");
    debugBuffer += L"\nВремя расчета: " + to_wstring(end) + L" миллисекунд\n";
    debugInfo().Text(debugBuffer);

    x_scale = canvas().ActualHeight() / Solver().H();
    y_scale = canvas().ActualWidth() / Solver().L();

    double xs = x_scale;
    double ys = y_scale;

    dx = Solver().Dx();
    dy = Solver().Dy();

    n = 0;

    for (int i = 0; i < Solver().Nx() * Solver().Ny(); i++)
    {
        double x = SolverImp()->nodes[i]->x;
        double y = SolverImp()->nodes[i]->y;
        //double maxV = SolverImp()->MaxVelocity(n);

        auto rect = Rect((x - dx / 2) * x_scale, canvas().ActualHeight() - (y + dy / 2) * y_scale, dx * x_scale, dy * y_scale);
        cells.push_back(rect);
    }
    
    canvas().Paused(false);
    
}

Windows::UI::Color winrt::CFD::implementation::CFDPage::ValueColor(double value)
{
    Color color = ColorHelper::FromArgb(255,
        (uint8_t)((value > 0.5 ? 2 * value - 1 : 0) * 255),
        (uint8_t)((value > 0.5 ? 2 - 2 * value : 2 * value) * 255),
        (uint8_t)((value > 0.5 ? 0 : (1 - 2 * value)) * 255));
    return color;
}


void winrt::CFD::implementation::CFDPage::canvas_Draw(winrt::Microsoft::Graphics::Canvas::UI::Xaml::ICanvasAnimatedControl const& sender, winrt::Microsoft::Graphics::Canvas::UI::Xaml::CanvasAnimatedDrawEventArgs const& args)
{
    double maxV = SolverImp()->MaxVelocity(n);
    for (int i = 0; i < Solver().Nx() * Solver().Ny(); i++)
    {
        Color col = ValueColor(SolverImp()->nodes[i]->AbsVelocity(n) / maxV);
        args.DrawingSession().FillRectangle(cells[i], col);
        //args.DrawingSession().FillEllipse(155, 115, 80, 30, col);
        //args.DrawingSession().DrawText("Hello, world!", 100, 100, Colors::Yellow);
    }
    n++;
    if (n == Solver().Nt() - 1)
        n = 0;
}


void winrt::CFD::implementation::CFDPage::Page_Loaded(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::RoutedEventArgs const& e)
{
    canvas().Paused(true);
}


void winrt::CFD::implementation::CFDPage::stopButton_Click(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::RoutedEventArgs const& e)
{
    solution.Cancel();
}


void winrt::CFD::implementation::CFDPage::Page_Unloaded(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::RoutedEventArgs const& e)
{
    canvas().RemoveFromVisualTree();
    canvas() = nullptr;
}
