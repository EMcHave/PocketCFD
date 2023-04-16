#pragma once

#include "CFDPage.g.h"
#include "NavierStokes.h"
#include <ppltasks.h>
#include <map>

using namespace winrt;
using namespace std;
using namespace Microsoft::Graphics::Canvas;
using namespace Microsoft::Graphics::Canvas::Effects;
using namespace Microsoft::Graphics::Canvas::Text;
using namespace Microsoft::Graphics::Canvas::UI::Xaml;
using namespace Microsoft::Graphics::Canvas::UI;
using namespace Windows::Foundation;
using namespace Windows::Foundation::Numerics;
using namespace Windows::Graphics::Effects;
using namespace Windows::Graphics::Imaging;
using namespace Windows::Storage;
using namespace Windows::Storage::Search;
using namespace Windows::Storage::Streams;
using namespace Windows::Storage::Pickers;
using namespace Windows::UI;
using namespace Windows::UI::Composition;
using namespace Windows::UI::Xaml;
using namespace Windows::UI::Xaml::Controls;
using namespace Windows::UI::Xaml::Input;
using namespace Windows::UI::Xaml::Media::Imaging;
using namespace Windows::UI::Xaml::Media::Animation;
using namespace Windows::UI::Xaml::Navigation;

namespace winrt::CFD::implementation
{
    struct CFDPage : CFDPageT<CFDPage>
    {
        CFDPage() 
        {
            Solver(make<NavierStokes>());
            Solver().Solved(false);
        }
        
        CFD::NavierStokes Solver() { return solver; }

        void Solver(CFD::NavierStokes const& value) { solver = value; }
        com_ptr<CFD::implementation::NavierStokes> SolverImp() { return Solver().as<CFD::implementation::NavierStokes>(); }

        IAsyncAction solveButton_Click(IInspectable const& sender, RoutedEventArgs const& e);

        Color ValueColor(double value);

        
    
    private:
        CFD::NavierStokes solver{ nullptr };
        wstring debugBuffer;
        vector<Rect> cells;
        vector<Rect> walls;
        vector<pair<float2, float2>> vectors;
        event_token m_showResidual;
        std::map<Boundary, BoundaryCondition> BoundaryConditions;
        vector<Point> boundaries;
        vector<Point> propants;
        BoundaryCondition left_bc, right_bc, top_bc, bottom_bc;

        Microsoft::Graphics::Canvas::Geometry::CanvasGeometry residualsPath_U{ nullptr };
        Microsoft::Graphics::Canvas::Geometry::CanvasGeometry residualsPath_V{ nullptr };
        Microsoft::Graphics::Canvas::Geometry::CanvasGeometry residualsPath_P{ nullptr };
        vector<double> residualsU, residualsV, residualsP;
        bool solving;
        bool animating;
        double Min;

        int n;
        int field;
        double dx;
        double dy;
        double x_scale;
        double y_scale;
        double height;
        double width;
        double dr;
        double resScale;
        IAsyncActionWithProgress<IVector<double>> solution{ nullptr };
    public:
        void ReDrawMesh();
        void PrepareAnimation();
        void PrepareResidualsPlot(Collections::IVector<double>& progress);
        void canvas_Draw(ICanvasAnimatedControl const& sender, CanvasAnimatedDrawEventArgs const& args);
        void Page_Loaded(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::RoutedEventArgs const& e);
        void stopButton_Click(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::RoutedEventArgs const& e);
        void Page_Unloaded(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::RoutedEventArgs const& e);
        fire_and_forget nxField_ValueChanged(winrt::Microsoft::UI::Xaml::Controls::NumberBox const& sender, winrt::Microsoft::UI::Xaml::Controls::NumberBoxValueChangedEventArgs const& args);
        fire_and_forget nyField_ValueChanged(winrt::Microsoft::UI::Xaml::Controls::NumberBox const& sender, winrt::Microsoft::UI::Xaml::Controls::NumberBoxValueChangedEventArgs const& args);
    
        void canvas_PointerPressed(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::Input::PointerRoutedEventArgs const& e);
        void Page_SizeChanged(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::SizeChangedEventArgs const& e);
        fire_and_forget lField_ValueChanged(winrt::Microsoft::UI::Xaml::Controls::NumberBox const& sender, winrt::Microsoft::UI::Xaml::Controls::NumberBoxValueChangedEventArgs const& args);
        fire_and_forget hField_ValueChanged(winrt::Microsoft::UI::Xaml::Controls::NumberBox const& sender, winrt::Microsoft::UI::Xaml::Controls::NumberBoxValueChangedEventArgs const& args);
        void leftBC_SelectionChanged(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::Controls::SelectionChangedEventArgs const& e);
        void rightBC_SelectionChanged(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::Controls::SelectionChangedEventArgs const& e);
        void topBC_SelectionChanged(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::Controls::SelectionChangedEventArgs const& e);
        void bottomBC_SelectionChanged(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::Controls::SelectionChangedEventArgs const& e);
        void leftSpeedValue_ValueChanged(winrt::Microsoft::UI::Xaml::Controls::NumberBox const& sender, winrt::Microsoft::UI::Xaml::Controls::NumberBoxValueChangedEventArgs const& args);
        void rightSpeedValue_ValueChanged(winrt::Microsoft::UI::Xaml::Controls::NumberBox const& sender, winrt::Microsoft::UI::Xaml::Controls::NumberBoxValueChangedEventArgs const& args);
        void topSpeedValue_ValueChanged(winrt::Microsoft::UI::Xaml::Controls::NumberBox const& sender, winrt::Microsoft::UI::Xaml::Controls::NumberBoxValueChangedEventArgs const& args);
        void bottomSpeedValue_ValueChanged(winrt::Microsoft::UI::Xaml::Controls::NumberBox const& sender, winrt::Microsoft::UI::Xaml::Controls::NumberBoxValueChangedEventArgs const& args);
        void fieldComboBox_SelectionChanged(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::Controls::SelectionChangedEventArgs const& e);
        fire_and_forget addBoundaryButton_Click(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::RoutedEventArgs const& e);
        fire_and_forget cleanBoundariesButton_Click(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::RoutedEventArgs const& e);
        void static_canvas_Draw(winrt::Microsoft::Graphics::Canvas::UI::Xaml::CanvasControl const& sender, winrt::Microsoft::Graphics::Canvas::UI::Xaml::CanvasDrawEventArgs const& args);
        void static_canvas_PointerMoved(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::Input::PointerRoutedEventArgs const& e);
        void animationButton_Click(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::RoutedEventArgs const& e);

        IAsyncOperation<ContentDialogResult> warningWindow()
        {
            ContentDialog dialog;
            dialog.Title(box_value(L"Внимание!"));
            dialog.Content(box_value(L"Изменение данного параметра приведет к очистке существующего решения"));
            dialog.PrimaryButtonText(L"OK");
            dialog.DefaultButton(ContentDialogButton::Primary);
            return dialog.ShowAsync();
        }
    };
}

namespace winrt::CFD::factory_implementation
{
    struct CFDPage : CFDPageT<CFDPage, implementation::CFDPage>
    {
    };
}
