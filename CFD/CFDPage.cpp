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
using namespace Windows::Foundation;

IAsyncAction winrt::CFD::implementation::CFDPage::solveButton_Click(IInspectable const& sender, RoutedEventArgs const& e)
{
    using namespace concurrency;

    canvas().Paused(true);

    progressBar().Visibility(Visibility::Visible());
    debugBuffer += L"Начат расчет\n";
    debugInfo().Text(debugBuffer);
    percentOfSolution().Text(L"0 %");

    int begin = GetTickCount();
    double T = Solver().T();
    solution = Solver().Solve();
    solution.Progress([&](auto const& sender, Collections::IVector<double> progress)
        {
            progressBar().Value(100*progress.GetAt(0));
            percentOfSolution().Text(to_wstring((int)(100*progress.GetAt(0))) + L" %");
            residualField().Value(progress.GetAt(1));
        });
    co_await solution;
    int end = GetTickCount() - begin;
    
    
    percentOfSolution().Text(L"Расчет завершен успешно");
    debugBuffer += L"\nВремя расчета: " + to_wstring(end) + L" миллисекунд\n";
    debugInfo().Text(debugBuffer);



    if (!cells.empty())
        cells.clear();

    x_scale = canvas().ActualWidth() / Solver().L();
    y_scale = canvas().ActualHeight() / Solver().H();

    double xs = x_scale;
    double ys = y_scale;

    dx = Solver().Dx();
    dy = Solver().Dy();

    float ratio = Solver().L() / Solver().H();
    if (ratio > 1)
    {
        canvas().Width(canvasBorder().ActualWidth() * 0.75);
        canvas().Height(canvas().Width() / ratio);
    }
    else
    {
        canvas().Height(canvasBorder().ActualHeight() * 0.75);
        canvas().Width(canvas().Height() * ratio);
    }

    for (int i = 0; i < Solver().Nx() * Solver().Ny(); i++)
    {
        double x = SolverImp()->nodes[i]->x;
        double y = SolverImp()->nodes[i]->y;

        auto rect = Rect((x - dx / 2) * x_scale, canvas().Height() - (y + dy / 2) * y_scale, dx * x_scale, dy * y_scale);
        cells.push_back(rect);
    }


    n = 0;
    Solver().Solved(true);
    canvas().Paused(false);
    
}

Color winrt::CFD::implementation::CFDPage::ValueColor(double value)
{
    Color color = ColorHelper::FromArgb(255,
        (uint8_t)((value > 0.5 ? 2 * value - 1 : 0) * 255),
        (uint8_t)((value > 0.5 ? 2 - 2 * value : 2 * value) * 255),
        (uint8_t)((value > 0.5 ? 0 : (1 - 2 * value)) * 255));
    return color;
}


void winrt::CFD::implementation::CFDPage::ReDrawMesh()
{
    if (!cells.empty())
        cells.clear();

    float ratio = Solver().L() / Solver().H();
    if (ratio > 1)
    {
        canvas().Width(canvasBorder().ActualWidth() * 0.75);
        canvas().Height(canvas().Width() / ratio);
    }   
    else
    {
        canvas().Height(canvasBorder().ActualHeight() * 0.75);
        canvas().Width(canvas().Height() * ratio);
    }
        

    double x, y;

    x_scale = canvas().Width() / Solver().L();
    y_scale = canvas().Height() / Solver().H();

    double xs = x_scale;
    double ys = y_scale;

    double dx = Solver().Dx();
    double dy = Solver().Dy();
    for (int i = 0; i < Solver().Nx(); i++)
        for (int j = 0; j < Solver().Ny(); j++)
        {
            x = i * dx;
            y = j * dy;
            auto rect = Rect(x * x_scale, canvas().Height() - y * y_scale, dx * x_scale, dy * y_scale);
            cells.push_back(rect);
        }
}

void winrt::CFD::implementation::CFDPage::canvas_Draw(ICanvasAnimatedControl const& sender, CanvasAnimatedDrawEventArgs const& args)
{
    if (Solver().Solved())
    {
        double maxV = SolverImp()->MaxVelocity(n);
        for (int i = 0; i < Solver().Nx() * Solver().Ny(); i++)
        {
            Color col = ValueColor(SolverImp()->nodes[i]->AbsVelocity(n) / maxV);
            args.DrawingSession().FillRectangle(cells[i], col);
        }
        n++;
        if (n == Solver().Nt() - 1)
            n = 0;
    }
    else
    {
        for (int i = 0; i < Solver().Nx() * Solver().Ny(); i++)
            args.DrawingSession().DrawRectangle(cells[i], Colors::Green());
    }
    
}


void winrt::CFD::implementation::CFDPage::Page_Loaded(IInspectable const& sender, RoutedEventArgs const& e)
{
    ReDrawMesh();
    canvas().Invalidate();
}


void winrt::CFD::implementation::CFDPage::stopButton_Click(IInspectable const& sender, RoutedEventArgs const& e)
{
    solution.Cancel();
}


void winrt::CFD::implementation::CFDPage::Page_Unloaded(IInspectable const& sender, RoutedEventArgs const& e)
{
    canvas().RemoveFromVisualTree();
    canvas() = nullptr;
}


void winrt::CFD::implementation::CFDPage::nxField_ValueChanged(winrt::Microsoft::UI::Xaml::Controls::NumberBox const& sender, winrt::Microsoft::UI::Xaml::Controls::NumberBoxValueChangedEventArgs const& args)
{
    Solver().Nx((int)nxField().Value());
    ReDrawMesh();
    canvas().Invalidate();
}


void winrt::CFD::implementation::CFDPage::nyField_ValueChanged(winrt::Microsoft::UI::Xaml::Controls::NumberBox const& sender, winrt::Microsoft::UI::Xaml::Controls::NumberBoxValueChangedEventArgs const& args)
{
    Solver().Ny((int)nyField().Value());
    ReDrawMesh();
    canvas().Invalidate();
}


void winrt::CFD::implementation::CFDPage::canvas_PointerPressed(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::Input::PointerRoutedEventArgs const& e)
{
    auto position = e.GetCurrentPoint(canvas()).Position();

    float i = std::round((position.X) / (Solver().Dx() * x_scale));
    float j = std::round((position.Y - canvas().ActualHeight()) / (-Solver().Dy() * y_scale));

    //debugBuffer += L"Скорость #" + to_wstring(SolverImp()->Nod(i, j)->u[n]) + L", " + to_wstring(SolverImp()->Nod(i, j)->v[n]) + L"\n";
    //debugBuffer += L"Выбран узел #" + to_wstring(SolverImp()->Nod(i, j)->ID) + L"\n";
    debugInfo().Text(L"Скорость #" + to_wstring(SolverImp()->Nod(i, j)->u[n]) + L", " + to_wstring(SolverImp()->Nod(i, j)->v[n]) + L"\n");

}


void winrt::CFD::implementation::CFDPage::Page_SizeChanged(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::SizeChangedEventArgs const& e)
{
    //ReDrawMesh();
    //canvas().Invalidate();
}


void winrt::CFD::implementation::CFDPage::lField_ValueChanged(winrt::Microsoft::UI::Xaml::Controls::NumberBox const& sender, winrt::Microsoft::UI::Xaml::Controls::NumberBoxValueChangedEventArgs const& args)
{
    Solver().L((double)lField().Value());
    ReDrawMesh();
    canvas().Invalidate();
}


void winrt::CFD::implementation::CFDPage::hField_ValueChanged(winrt::Microsoft::UI::Xaml::Controls::NumberBox const& sender, winrt::Microsoft::UI::Xaml::Controls::NumberBoxValueChangedEventArgs const& args)
{
    Solver().H((double)hField().Value());
    ReDrawMesh();
    canvas().Invalidate();
}
