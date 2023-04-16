#include "pch.h"
#include "CFDPage.h"
#include <chrono>
#include <Windows.h>
#include <unknwn.h>
#include <pplawait.h>
#if __has_include("CFDPage.g.cpp")
#include "CFDPage.g.cpp"
#endif

using namespace winrt;
using namespace Windows::UI;
using namespace Windows::UI::Xaml;
using namespace winrt::Microsoft::UI::Xaml::Controls;
using namespace Microsoft::Graphics::Canvas::UI::Xaml;
using namespace Windows::Foundation;

IAsyncAction winrt::CFD::implementation::CFDPage::solveButton_Click(IInspectable const& sender, RoutedEventArgs const& e)
{
    using namespace concurrency;

    canvas().Paused(false);
    Solver().Solved(false);
    ReDrawMesh();
    canvas().Visibility(Visibility::Visible);

    
    progressBar().Visibility(Visibility::Visible());
    debugBuffer += L"Начат расчет\n";
    debugInfo().Text(debugBuffer);
    percentOfSolution().Text(L"0 %");

    int begin = GetTickCount64();

    BoundaryConditions[Boundary::Left] = left_bc;
    BoundaryConditions[Boundary::Right] = right_bc;
    BoundaryConditions[Boundary::Top] = top_bc;
    BoundaryConditions[Boundary::Bottom] = bottom_bc;

    solving = true;

    residualsU.clear();
    residualsV.clear();
    residualsP.clear();
    
    int itCount = 0;
    
    Min = Solver().Eps();

    solution = SolverImp()->Solve(BoundaryConditions, boundaries);
    solution.Progress([&](auto const& sender, Collections::IVector<double> progress)
        {
            progressBar().Value(100 * progress.GetAt(0));
            percentOfSolution().Text(to_wstring((int)(100 * progress.GetAt(0))) + L" %");
            residualField().Value(progress.GetAt(4));

            if (itCount > 0)
                PrepareResidualsPlot(progress);
            itCount++;
        });
    co_await solution;
    int end = GetTickCount64() - begin;
    
    solving = false;
    
    percentOfSolution().Text(L"Расчет завершен успешно");
    debugBuffer += L"\nВремя расчета: " + to_wstring(end) + L" миллисекунд\n";
    debugInfo().Text(debugBuffer);
    

    
    Solver().Solved(true);
    n = 0;
    PrepareAnimation();
    canvas().Paused(false);
}

void winrt::CFD::implementation::CFDPage::ReDrawMesh()
{
    if (!cells.empty())
        cells.clear();

    float ratio = Solver().L() / Solver().H();
    if (ratio > 1)
    {
        static_canvas().Width(canvasBorder().ActualWidth() * 0.75);
        static_canvas().Height(static_canvas().Width() / ratio);
    }
    else
    {
        static_canvas().Height(canvasBorder().ActualHeight() * 0.75);
        static_canvas().Width(static_canvas().Height() * ratio);
    }

    canvas().Visibility(Visibility::Collapsed);
    static_canvas().Visibility(Visibility::Visible);

    height = static_canvas().Height();
    width = static_canvas().Width();

    x_scale = static_canvas().Width() / Solver().L();
    y_scale = static_canvas().Height() / Solver().H();

    double dx = Solver().Dx();
    double dy = Solver().Dy();

    double x, y;
    for (int i = 0; i < Solver().Nx(); i++)
        for (int j = 0; j < Solver().Ny(); j++)
        {
            x = i * dx;
            y = j * dy;
            auto rect = Rect(x * x_scale, static_canvas().Height() - y * y_scale, dx * x_scale, dy * y_scale);
            cells.push_back(rect);
        }

    for (Point p : boundaries)
    {
        double dx = Solver().Dx();
        double dy = Solver().Dy();
        auto rect = Rect(p.X*dx * x_scale, static_canvas().Height() - p.Y*dy * y_scale, dx * x_scale, dy * y_scale);
        walls.push_back(rect);
    }
}

void winrt::CFD::implementation::CFDPage::PrepareAnimation()
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

    canvas().Visibility(Visibility::Visible);
    static_canvas().Visibility(Visibility::Collapsed);

    height = canvas().Height();
    width = canvas().Width();

    x_scale = canvas().ActualWidth() / Solver().L();
    y_scale = canvas().ActualHeight() / Solver().H();

    dx = Solver().Dx();
    dy = Solver().Dy(); 

    switch (field)
    {
    case 0:
    {
        for (int i = 0; i < Solver().Nx() * Solver().Ny(); i++)
        {
            double x = SolverImp()->nodes[i]->x;
            double y = SolverImp()->nodes[i]->y;

            auto rect = Rect((x - dx / 2) * x_scale, (y - dy / 2) * y_scale, dx * x_scale, dy * y_scale);
            cells.push_back(rect);
        }
        break;
    }
    case 1:
    {
        for (int i = 0; i < (Solver().Nx() - 1) * (Solver().Ny() - 1); i++)
        {

            double x = SolverImp()->cells[i]->n1->x;
            double y = SolverImp()->cells[i]->n4->y;

            auto rect = Rect(x * x_scale, y * y_scale, dx * x_scale, dy * y_scale);
            cells.push_back(rect);
        }
        break;
    }
    default:
        break;
    }
}

void winrt::CFD::implementation::CFDPage::PrepareResidualsPlot(Collections::IVector<double>& progress)
{
    double V0 = progress.GetAt(2);
    residualsU.push_back(progress.GetAt(1));
    residualsV.push_back(progress.GetAt(2));
    residualsP.push_back(progress.GetAt(3));
    double max = max(residualsU[0], max(residualsV[0], residualsP[0]));


    Min = min(Min, min(residualsU.back(), min(residualsV.back(), residualsP.back())));

    dr = log10(max) - log10(Min / 10);
    resScale = static_canvas().Height() / dr;

    Microsoft::Graphics::Canvas::Geometry::CanvasPathBuilder residualBuild_U(static_canvas().as<ICanvasResourceCreator>());
    Microsoft::Graphics::Canvas::Geometry::CanvasPathBuilder residualBuild_V(static_canvas().as<ICanvasResourceCreator>());
    Microsoft::Graphics::Canvas::Geometry::CanvasPathBuilder residualBuild_P(static_canvas().as<ICanvasResourceCreator>());


    residualBuild_U.BeginFigure(0, static_canvas().Height() - resScale * (log10(residualsU[0]) - log10(Min / 10)));
    residualBuild_V.BeginFigure(0, static_canvas().Height() - resScale * (log10(residualsV[0]) - log10(Min / 10)));
    residualBuild_P.BeginFigure(0, static_canvas().Height() - resScale * (log10(residualsP[0]) - log10(Min / 10)));

    for (int it = 1; it < residualsU.size(); it++)
    {
        double xiu = static_canvas().Height() - resScale * (log10(residualsU[it]) - log10(Min / 10));
        double xiv = static_canvas().Height() - resScale * (log10(residualsV[it]) - log10(Min / 10));
        residualBuild_U.AddLine(it * static_canvas().Width() / residualsU.size(), xiu);
        residualBuild_V.AddLine(it * static_canvas().Width() / residualsU.size(), xiv);
        residualBuild_P.AddLine(it * static_canvas().Width() / residualsP.size(), static_canvas().Height() - resScale * (log10(residualsP[it]) - log10(Min / 10)));
    }

    residualBuild_U.EndFigure(Microsoft::Graphics::Canvas::Geometry::CanvasFigureLoop::Open);
    residualBuild_V.EndFigure(Microsoft::Graphics::Canvas::Geometry::CanvasFigureLoop::Open);
    residualBuild_P.EndFigure(Microsoft::Graphics::Canvas::Geometry::CanvasFigureLoop::Open);
    residualsPath_U = Microsoft::Graphics::Canvas::Geometry::CanvasGeometry::CreatePath(residualBuild_U);
    residualsPath_V = Microsoft::Graphics::Canvas::Geometry::CanvasGeometry::CreatePath(residualBuild_V);
    residualsPath_P = Microsoft::Graphics::Canvas::Geometry::CanvasGeometry::CreatePath(residualBuild_P);


    static_canvas().Invalidate();
}

Color winrt::CFD::implementation::CFDPage::ValueColor(double value)
{
    Color color = ColorHelper::FromArgb(255,
        (uint8_t)((value > 0.5 ? 2 * value - 1 : 0) * 255),
        (uint8_t)((value > 0.5 ? 2 - 2 * value : 2 * value) * 255),
        (uint8_t)((value > 0.5 ? 0 : (1 - 2 * value)) * 255));
    return color;
}

void winrt::CFD::implementation::CFDPage::canvas_Draw(ICanvasAnimatedControl const& sender, CanvasAnimatedDrawEventArgs const& args)
{
    if (Solver().Solved())
    {
        double maxV = 0.0;
        double maxP, minP;
        int N = 0;
        switch (field)
        {
        case 0:
        {
            maxV = SolverImp()->MaxVelocity(n);
            N = Solver().Nx() * Solver().Ny();
            for (int i = 0; i < N; i++)
            {
                Color col = ValueColor(SolverImp()->nodes[i]->AbsVelocity(n) / maxV);
                args.DrawingSession().FillRectangle(cells[i], col);
            }
            break;
        }
        case 1:
        {
            maxP = SolverImp()->MaxPressure(n);
            minP = SolverImp()->MinPressure(n);
            N = (Solver().Nx() - 1) * (Solver().Ny() - 1);
            for (int i = 0; i < N; i++)
            {
                Color col = ValueColor( (SolverImp()->cells[i]->p[n] - minP) / (maxP - minP));
                args.DrawingSession().FillRectangle(cells[i], col);
            }
            break;
        }
        case 2:
        {
            maxV = SolverImp()->MaxVelocity(n);
            N = Solver().Nx() * Solver().Ny();
            for (int i = 0; i < N; i++)
            {
                float2 f1;
                float2 f2;
                Color col = ValueColor(SolverImp()->nodes[i]->AbsVelocity(n) / maxV);
                double x = SolverImp()->nodes[i]->x;
                double y = SolverImp()->nodes[i]->y;
                double u = SolverImp()->nodes[i]->u[n];
                double v = SolverImp()->nodes[i]->v[n];
                f1.x = x * x_scale;
                f1.y = y * y_scale;
                f2.x = (x + u / maxV * SolverImp()->Dx()) * x_scale;
                f2.y = (y - v / maxV * SolverImp()->Dy()) * y_scale;
                args.DrawingSession().DrawLine(f1, f2, col, 2);
            }
            break;
        }
        default:
            break;
        }

        n++;
        if (n == Solver().Nt() - 1)
            n = 0;
    }
    if (animating)
    {
        for (int i = 0; i < propants.size(); i++)
        {
            Point p = propants[i];
            int nod_i = (int)(propants[i].X / Solver().Dx());
            int nod_j = (int)(propants[i].Y / Solver().Dy());
            double x_new = propants[i].X + SolverImp()->Nod((int)(propants[i].X / Solver().Dx()), (int)(propants[i].Y / Solver().Dy()))->u[n] * Solver().Dt();
            double y_new = propants[i].Y + SolverImp()->Nod((int)(propants[i].X / Solver().Dx()), (int)(propants[i].Y / Solver().Dy()))->v[n] * Solver().Dt();

            int i_new = x_new / Solver().Dx();
            int j_new = y_new / Solver().Dy();

            if (i_new <= Solver().Nx() && j_new <= Solver().Ny())
            {
                propants[i].X = x_new;
                propants[i].Y = y_new;
                args.DrawingSession().FillCircle(x_new * x_scale, height - y_new * y_scale, max(Solver().Dy() / 2, Solver().Dx() / 2) * max(x_scale, y_scale) * 0.9, Colors::Coral());
            }
            else
                propants.erase(propants.begin() + i);
        }
        if (propants.empty())
            animating = false;
    }
    
}

void winrt::CFD::implementation::CFDPage::fieldComboBox_SelectionChanged(IInspectable const& sender, SelectionChangedEventArgs const& e)
{
    canvas().Paused(true);
    field = sender.as<ComboBox>().SelectedIndex();
    
    if (Solver().Solved())
    {
        n = 0;
        PrepareAnimation();
        canvas().Paused(false);
    }
}

void winrt::CFD::implementation::CFDPage::Page_Loaded(IInspectable const& sender, RoutedEventArgs const& e)
{

    ReDrawMesh();
    static_canvas().Invalidate();
    leftBC().SelectedIndex(1);
    rightBC().SelectedIndex(2);
    topBC().SelectedIndex(0);
    bottomBC().SelectedIndex(0);
    fieldComboBox().SelectedIndex(0);
    walls = vector<Rect>();
    boundaries = vector<Point>();
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


fire_and_forget winrt::CFD::implementation::CFDPage::nxField_ValueChanged(NumberBox const& sender, NumberBoxValueChangedEventArgs const& args)
{
    if (Solver().Solved() == true)
    {
        auto result = warningWindow();
        Solver().Solved(false);
    }

    Solver().Nx((int)nxField().Value());
    ReDrawMesh();
    static_canvas().Invalidate();
}


fire_and_forget winrt::CFD::implementation::CFDPage::nyField_ValueChanged(NumberBox const& sender, NumberBoxValueChangedEventArgs const& args)
{
    if (Solver().Solved() == true)
    {
        auto result = warningWindow();
        Solver().Solved(false);
    }

    Solver().Ny((int)nyField().Value());
    ReDrawMesh();
    static_canvas().Invalidate();
}


void winrt::CFD::implementation::CFDPage::canvas_PointerPressed(IInspectable const& sender, winrt::Windows::UI::Xaml::Input::PointerRoutedEventArgs const& e)
{
    auto position = e.GetCurrentPoint(canvas()).Position();

    float i = std::round((position.X) / (Solver().Dx() * x_scale));
    float j = std::round((position.Y - canvas().ActualHeight()) / (-Solver().Dy() * y_scale));
    if(Solver().Solved())
        switch (field)
        {
        case 0:
            if (i < Solver().Nx() && j < Solver().Ny())
                debugInfo().Text(L"Скорость #" + to_wstring(SolverImp()->Nod(i, j)->u[n]) + L", " + to_wstring(SolverImp()->Nod(i, j)->v[n]) + L"\n");
            break;
        case 1:
            if(i < Solver().Nx() - 1 && j < Solver().Ny() - 1)
                debugInfo().Text(L"Давление #" + to_wstring(SolverImp()->Cel(i, j)->p[n]) + L"\n");
            break;
        default:
            break;
        }
}


void winrt::CFD::implementation::CFDPage::Page_SizeChanged(IInspectable const& sender, winrt::Windows::UI::Xaml::SizeChangedEventArgs const& e)
{
    //ReDrawMesh();
    //canvas().Invalidate();
}


fire_and_forget winrt::CFD::implementation::CFDPage::lField_ValueChanged(NumberBox const& sender, NumberBoxValueChangedEventArgs const& args)
{
    if (Solver().Solved() == true)
    {
        auto result = warningWindow();
        Solver().Solved(false);
    }

    Solver().L((double)lField().Value());
    ReDrawMesh();
    static_canvas().Invalidate();
}


fire_and_forget winrt::CFD::implementation::CFDPage::hField_ValueChanged(NumberBox const& sender, NumberBoxValueChangedEventArgs const& args)
{
    if (Solver().Solved() == true)
    {
        auto result = warningWindow();
        Solver().Solved(false);
    }

    Solver().H((double)hField().Value());
    ReDrawMesh();
    static_canvas().Invalidate();
}


void winrt::CFD::implementation::CFDPage::leftBC_SelectionChanged(IInspectable const& sender, SelectionChangedEventArgs const& e)
{
    if (leftBC().SelectedIndex() == 0)
    {
        leftSpeed().Visibility(Visibility::Visible);
        leftPressure().Visibility(Visibility::Collapsed);
        leftSpeedXValue().Value(0);
        leftSpeedXValue().IsEnabled(false);
        leftSpeedYValue().Value(0);
        leftSpeedYValue().IsEnabled(false);
    }
    else if (leftBC().SelectedIndex() == 1)
    {
        leftSpeed().Visibility(Visibility::Visible);
        leftPressure().Visibility(Visibility::Collapsed);
        leftSpeedXValue().Value(1.0);
        leftSpeedXValue().IsEnabled(true);
        leftSpeedYValue().Value(0);
        leftSpeedYValue().IsEnabled(true);
    }
    else if (leftBC().SelectedIndex() == 2)
    {
        leftSpeed().Visibility(Visibility::Collapsed);
        leftPressure().Visibility(Visibility::Visible);
        leftSpeedXValue().Value(0);
        leftSpeedXValue().IsEnabled(false);
        leftSpeedYValue().Value(0);
        leftSpeedYValue().IsEnabled(false);
    }
    left_bc.Type = static_cast<TypeOfBC>(sender.as<ComboBox>().SelectedIndex());
    leftPressureValue().Value(0);
}


void winrt::CFD::implementation::CFDPage::rightBC_SelectionChanged(IInspectable const& sender, SelectionChangedEventArgs const& e)
{
    if (rightBC().SelectedIndex() == 0)
    {
        rightSpeed().Visibility(Visibility::Visible);
        rightPressure().Visibility(Visibility::Collapsed);
        rightSpeedXValue().Value(0);
        rightSpeedXValue().IsEnabled(false);
        rightSpeedYValue().Value(0);
        rightSpeedYValue().IsEnabled(false);
    }
    else if (rightBC().SelectedIndex() == 1)
    {
        rightSpeed().Visibility(Visibility::Visible);
        rightPressure().Visibility(Visibility::Collapsed);
        rightSpeedXValue().Value(1.0);
        rightSpeedXValue().IsEnabled(true);
        rightSpeedYValue().Value(0);
        rightSpeedYValue().IsEnabled(true);
    }
    else if (rightBC().SelectedIndex() == 2)
    {
        rightSpeed().Visibility(Visibility::Collapsed);
        rightPressure().Visibility(Visibility::Visible);
        rightSpeedXValue().Value(0);
        rightSpeedXValue().IsEnabled(false);
        rightSpeedYValue().Value(0);
        rightSpeedYValue().IsEnabled(false);
    }
    right_bc.Type = static_cast<TypeOfBC>(sender.as<ComboBox>().SelectedIndex());
    rightPressureValue().Value(0);
}


void winrt::CFD::implementation::CFDPage::topBC_SelectionChanged(IInspectable const& sender, SelectionChangedEventArgs const& e)
{
    if (topBC().SelectedIndex() == 0)
    {
        topSpeed().Visibility(Visibility::Visible);
        topPressure().Visibility(Visibility::Collapsed);
        topSpeedXValue().Value(0);
        topSpeedXValue().IsEnabled(false);
        topSpeedYValue().Value(0);
        topSpeedYValue().IsEnabled(false);
    }
    else if (topBC().SelectedIndex() == 1)
    {
        topSpeed().Visibility(Visibility::Visible);
        topPressure().Visibility(Visibility::Collapsed);
        topSpeedXValue().Value(0);
        topSpeedXValue().IsEnabled(true);
        topSpeedYValue().Value(1.0);
        topSpeedYValue().IsEnabled(true);
    }
    else if (topBC().SelectedIndex() == 2)
    {
        topSpeed().Visibility(Visibility::Collapsed);
        topPressure().Visibility(Visibility::Visible);
        topSpeedXValue().Value(0);
        topSpeedXValue().IsEnabled(false);
        topSpeedYValue().Value(0);
        topSpeedYValue().IsEnabled(false);
    }
    top_bc.Type = static_cast<TypeOfBC>(sender.as<ComboBox>().SelectedIndex());
    topPressureValue().Value(0);
}


void winrt::CFD::implementation::CFDPage::bottomBC_SelectionChanged(IInspectable const& sender, SelectionChangedEventArgs const& e)
{
    if (bottomBC().SelectedIndex() == 0)
    {
        bottomSpeed().Visibility(Visibility::Visible);
        bottomPressure().Visibility(Visibility::Collapsed);
        bottomSpeedXValue().Value(0);
        bottomSpeedXValue().IsEnabled(false);
        bottomSpeedYValue().Value(0);
        bottomSpeedYValue().IsEnabled(false);
    }
    else if (bottomBC().SelectedIndex() == 1)
    {
        bottomSpeed().Visibility(Visibility::Visible);
        bottomPressure().Visibility(Visibility::Collapsed);
        bottomSpeedXValue().Value(0);
        bottomSpeedXValue().IsEnabled(true);
        bottomSpeedYValue().Value(1.0);
        bottomSpeedYValue().IsEnabled(true);
    }
    else if (bottomBC().SelectedIndex() == 2)
    {
        bottomSpeed().Visibility(Visibility::Collapsed);
        bottomPressure().Visibility(Visibility::Visible);
        bottomSpeedXValue().Value(0);
        bottomSpeedXValue().IsEnabled(false);
        bottomSpeedYValue().Value(0);
        bottomSpeedYValue().IsEnabled(false);
    }
    bottom_bc.Type = static_cast<TypeOfBC>(sender.as<ComboBox>().SelectedIndex());
    bottomPressureValue().Value(0);
}


void winrt::CFD::implementation::CFDPage::leftSpeedValue_ValueChanged(NumberBox const& sender, NumberBoxValueChangedEventArgs const& args)
{
    left_bc.Vx = leftSpeedXValue().Value();
    left_bc.Vy = leftSpeedYValue().Value();
    left_bc.P = leftPressureValue().Value();
}


void winrt::CFD::implementation::CFDPage::rightSpeedValue_ValueChanged(NumberBox const& sender, NumberBoxValueChangedEventArgs const& args)
{
    right_bc.Vx = rightSpeedXValue().Value();
    right_bc.Vy = rightSpeedYValue().Value();
    right_bc.P = rightPressureValue().Value();
}


void winrt::CFD::implementation::CFDPage::topSpeedValue_ValueChanged(NumberBox const& sender, NumberBoxValueChangedEventArgs const& args)
{
    top_bc.Vx = topSpeedXValue().Value();
    top_bc.Vy = topSpeedYValue().Value();
    top_bc.P = topPressureValue().Value();
}


void winrt::CFD::implementation::CFDPage::bottomSpeedValue_ValueChanged(NumberBox const& sender, NumberBoxValueChangedEventArgs const& args)
{
    bottom_bc.Vx = bottomSpeedXValue().Value();
    bottom_bc.Vy = bottomSpeedYValue().Value();
    bottom_bc.P = bottomPressureValue().Value();
}


fire_and_forget winrt::CFD::implementation::CFDPage::addBoundaryButton_Click(IInspectable const& sender, winrt::Windows::UI::Xaml::RoutedEventArgs const& e)
{


    if (true)
    {
        for (int i = (int)n1_i_Field().Value(); i <= (int)n2_i_Field().Value(); i++)
            for (int j = (int)n1_j_Field().Value(); j <= (int)n2_j_Field().Value(); j++)
                boundaries.push_back(Point(i, j));
    }
    ReDrawMesh();
    static_canvas().Invalidate();
}


fire_and_forget winrt::CFD::implementation::CFDPage::cleanBoundariesButton_Click(IInspectable const& sender, winrt::Windows::UI::Xaml::RoutedEventArgs const& e)
{


    boundaries.clear();
    walls.clear();
    ReDrawMesh();
    static_canvas().Invalidate();
}


void winrt::CFD::implementation::CFDPage::static_canvas_Draw(winrt::Microsoft::Graphics::Canvas::UI::Xaml::CanvasControl const& sender, winrt::Microsoft::Graphics::Canvas::UI::Xaml::CanvasDrawEventArgs const& args)
{
    if (solving)
    {
        args.DrawingSession().DrawGeometry(residualsPath_U, Colors::Red());
        args.DrawingSession().DrawGeometry(residualsPath_V, Colors::Green());
        args.DrawingSession().DrawGeometry(residualsPath_P, Colors::Cyan());
        args.DrawingSession().DrawText(L"X-Velocity", float2(width - 100, 10), Colors::Red());
        args.DrawingSession().DrawText(L"Y-Velocity", float2(width - 100, 40), Colors::Green());
        args.DrawingSession().DrawText(L"Continuity", float2(width - 100, 70), Colors::Cyan());
        args.DrawingSession().DrawLine(float2(0, height - resScale * (log10(Solver().Eps()) - log10(Min / 10))), float2(width, height - resScale * (log10(Solver().Eps()) - log10(Min / 10))), Colors::Orange());
        args.DrawingSession().DrawText(L"Eps", float2(0, height - resScale * (log10(Solver().Eps()) - log10(Min / 10)) - 30), Colors::Orange());
    }
    else
    {
        for (int i = 0; i < Solver().Nx() * Solver().Ny(); i++)
            args.DrawingSession().DrawRectangle(cells[i], Colors::Green());
        for (int i = 0; i < walls.size(); i++)
            args.DrawingSession().FillRectangle(walls[i], Colors::Green());
    }
}


void winrt::CFD::implementation::CFDPage::static_canvas_PointerMoved(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::Input::PointerRoutedEventArgs const& e)
{
    auto ptr = e.Pointer();

    float px = -1;
    float py = -1;


    if (ptr.PointerDeviceType() == Windows::Devices::Input::PointerDeviceType::Mouse)
    {
        Windows::UI::Input::PointerPoint ptrPt = e.GetCurrentPoint(static_canvas());
        if (ptrPt.Properties().IsLeftButtonPressed())
        {
            Point p;
            p.X = std::round((ptrPt.Position().X) / (Solver().Dx() * x_scale));
            p.Y = std::round((static_canvas().ActualHeight() - ptrPt.Position().Y) / (Solver().Dy() * y_scale));
            if (p.X != px && p.Y != py)
            {
                boundaries.push_back(p);
                px = p.X;
                py = p.Y;
                ReDrawMesh();
                static_canvas().Invalidate();
            }
        }

    }

}


void winrt::CFD::implementation::CFDPage::animationButton_Click(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::RoutedEventArgs const& e)
{
    if (Solver().Solved())
    {
        if (!propants.empty())
            propants.clear();

        switch (entryPointComboBox().SelectedIndex())
        {
        case 0:
            for (int i = 0; i < Solver().Ny(); i++)
                propants.push_back(Point(0, i * Solver().Dy()));
            break;
        case 1:
            for (int i = 0; i < Solver().Ny(); i++)
                propants.push_back(Point(Solver().Nx() * Solver().Dx(), i * Solver().Dy()));
            break;
        case 2:
            for (int i = 0; i < Solver().Nx(); i++)
                propants.push_back(Point(i * Solver().Dx(), Solver().Ny() * Solver().Dy()));
            break;
        case 3:
            for (int i = 0; i < Solver().Nx(); i++)
                propants.push_back(Point(i * Solver().Dx(), 0));
            break;

        default:
            break;
        }
        animating = true;
    }
}
