#pragma once

#include "CFDPage.g.h"
#include "NavierStokes.h"
#include <ppltasks.h>
#include <map>

#define MY_BLACK       float3(0.,0.,0.)   
#define MY_RED         float3(0.93, 0   , 0.01)
#define MY_LIGHT_RED   float3(1,  0.5 , 0.5)
#define MY_ORANGE      float3(1   , 0.45, 0.21)
#define MY_YELLOW      float3(0.9 , 0.9 , 0   )
#define MY_LIGHT_YELLOW float3(1 , 1 , 0.8 )  
#define MY_LIGHT_ROSE  float3(1 , 0.8 , 0.8 )  
#define MY_DARK_YELLOW float3(0.45 , 0.45, 0   )
#define MY_GREEN       float3(0.01, 0.98, 0.01)
#define MY_BLUE        float3(0   , 0.1 , 1   )
#define MY_LIGHT_BLUE  float3(0.5 , 0.75 , 1   )
#define MY_INDIGO      float3(0.3 , 0   , 0.53)
#define MY_VIOLET      float3(0.58, 0   , 0.83)
#define MY_DARK_RED    float3(0.5, 0   , 0.01)
#define MY_ROSE        float3(1., 0.  , 1.)
#define MY_PURPLE      float3(0.6, 0  , 0.6)
#define MY_DARK_GREEN  float3(0.01, 0.5, 0.01)
#define MY_DARK_BLUE   float3(0   , 0.1 , 0.5   )
#define MY_GOLD        float3(0.864 , 0.7 , 0.325 )
#define MY_DARK_GOLD   float3(0.4,    0.3 , 0.16 )
#define MY_SILVER      float3(0.75,0.75,0.75)
#define MY_PEATCH      float3(1, 0.89, 0.705)
#define MY_BRONZE      float3(0.8, 0.5, 0.2)
#define MY_GRAY        float3(0.5,0.5,0.5) 
#define MY_LIGHT_GRAY  float3(0.8,0.8,0.8) 
#define MY_DARK_GRAY   float3(0.2,0.2,0.2) 
#define MY_WHITE       float3(1.,1.,1.)
#define MY_DARK_CHERRY float3(0.57,0.11,0.26)
#define MY_LIGHT_CHERRY float3(0.87,0.2,0.4)



using namespace winrt;
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
        Color ValueColorBessonov(double x, int IsStriped, double scale);
        IAsyncAction ClearSolution();
    
    private:
        CFD::NavierStokes solver{ nullptr };
        std::wstring debugBuffer;
        std::vector<Rect> cells;
        std::vector<Rect> walls;
        std::vector<std::pair<float2, float2>> vectors;
        event_token m_showResidual;
        std::map<Boundary, BoundaryCondition> BoundaryConditions;
        std::vector<Point> boundaries;
        std::vector<Point> propants;
        BoundaryCondition left_bc, right_bc, top_bc, bottom_bc;
        const float3 RainBow7[7] = {MY_BLUE, MY_LIGHT_BLUE, MY_GREEN, MY_YELLOW, MY_ORANGE, MY_RED, MY_INDIGO};
        Windows::UI::Core::CoreDispatcher disp{ nullptr };

        Microsoft::Graphics::Canvas::Geometry::CanvasGeometry residualsPath_U{ nullptr };
        Microsoft::Graphics::Canvas::Geometry::CanvasGeometry residualsPath_V{ nullptr };
        Microsoft::Graphics::Canvas::Geometry::CanvasGeometry residualsPath_P{ nullptr };
        std::vector<double> residualsU, residualsV, residualsP;
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
        int sdvig;
        int colorScaleHeight;
        double dr;
        double resScale;
        double colorScale;
        IAsyncActionWithProgress<IVector<double>> solution{ nullptr };
    public:
        void ReDrawMesh();
        void PrepareAnimation();
        void PrepareResidualsPlot(Collections::IVector<double>& progress);
        void canvas_Draw(ICanvasAnimatedControl const& sender, CanvasAnimatedDrawEventArgs const& args);
        void Page_Loaded(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::RoutedEventArgs const& e);
        void stopButton_Click(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::RoutedEventArgs const& e);
        void Page_Unloaded(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::RoutedEventArgs const& e);
        IAsyncAction nxField_ValueChanged(winrt::Microsoft::UI::Xaml::Controls::NumberBox const& sender, winrt::Microsoft::UI::Xaml::Controls::NumberBoxValueChangedEventArgs const& args);
        IAsyncAction nyField_ValueChanged(winrt::Microsoft::UI::Xaml::Controls::NumberBox const& sender, winrt::Microsoft::UI::Xaml::Controls::NumberBoxValueChangedEventArgs const& args);
    
        void canvas_PointerPressed(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::Input::PointerRoutedEventArgs const& e);
        void Page_SizeChanged(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::SizeChangedEventArgs const& e);
        IAsyncAction lField_ValueChanged(winrt::Microsoft::UI::Xaml::Controls::NumberBox const& sender, winrt::Microsoft::UI::Xaml::Controls::NumberBoxValueChangedEventArgs const& args);
        IAsyncAction hField_ValueChanged(winrt::Microsoft::UI::Xaml::Controls::NumberBox const& sender, winrt::Microsoft::UI::Xaml::Controls::NumberBoxValueChangedEventArgs const& args);
        void leftBC_SelectionChanged(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::Controls::SelectionChangedEventArgs const& e);
        void rightBC_SelectionChanged(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::Controls::SelectionChangedEventArgs const& e);
        void topBC_SelectionChanged(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::Controls::SelectionChangedEventArgs const& e);
        void bottomBC_SelectionChanged(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::Controls::SelectionChangedEventArgs const& e);
        void leftSpeedValue_ValueChanged(winrt::Microsoft::UI::Xaml::Controls::NumberBox const& sender, winrt::Microsoft::UI::Xaml::Controls::NumberBoxValueChangedEventArgs const& args);
        void rightSpeedValue_ValueChanged(winrt::Microsoft::UI::Xaml::Controls::NumberBox const& sender, winrt::Microsoft::UI::Xaml::Controls::NumberBoxValueChangedEventArgs const& args);
        void topSpeedValue_ValueChanged(winrt::Microsoft::UI::Xaml::Controls::NumberBox const& sender, winrt::Microsoft::UI::Xaml::Controls::NumberBoxValueChangedEventArgs const& args);
        void bottomSpeedValue_ValueChanged(winrt::Microsoft::UI::Xaml::Controls::NumberBox const& sender, winrt::Microsoft::UI::Xaml::Controls::NumberBoxValueChangedEventArgs const& args);
        void fieldComboBox_SelectionChanged(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::Controls::SelectionChangedEventArgs const& e);
        IAsyncAction addBoundaryButton_Click(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::RoutedEventArgs const& e);
        IAsyncAction cleanBoundariesButton_Click(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::RoutedEventArgs const& e);
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
        void PauseButton_Checked(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::RoutedEventArgs const& e);
        void PauseButton_Unchecked(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::RoutedEventArgs const& e);
        void PauseButton_Click_1(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::RoutedEventArgs const& e);
        void colorSilder_ValueChanged(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::Controls::Primitives::RangeBaseValueChangedEventArgs const& e);
    };
}

namespace winrt::CFD::factory_implementation
{
    struct CFDPage : CFDPageT<CFDPage, implementation::CFDPage>
    {
    };
}
