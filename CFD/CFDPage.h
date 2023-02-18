#pragma once

#include "CFDPage.g.h"
#include "NavierStokes.h"
#include <ppltasks.h>

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
            //anvas().Paused(true);
        }
        
        CFD::NavierStokes Solver() { return solver; }

        void Solver(CFD::NavierStokes const& value) { solver = value; }
        com_ptr<CFD::implementation::NavierStokes> SolverImp() { return Solver().as<CFD::implementation::NavierStokes>(); }

        IAsyncAction solveButton_Click(IInspectable const& sender, RoutedEventArgs const& e);

        Color ValueColor(double value);


    
    private:
        CFD::NavierStokes solver{ nullptr };
        wstring debugBuffer;
        std::vector<Rect> cells;
        int n;
        double dx;
        double dy;
        double x_scale;
        double y_scale;
        IAsyncActionWithProgress<double> solution{ nullptr };
    public:
        void canvas_Draw(ICanvasAnimatedControl const& sender, CanvasAnimatedDrawEventArgs const& args);
        void Page_Loaded(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::RoutedEventArgs const& e);
        void stopButton_Click(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::RoutedEventArgs const& e);
        void Page_Unloaded(winrt::Windows::Foundation::IInspectable const& sender, winrt::Windows::UI::Xaml::RoutedEventArgs const& e);
    };
}

namespace winrt::CFD::factory_implementation
{
    struct CFDPage : CFDPageT<CFDPage, implementation::CFDPage>
    {
    };
}
