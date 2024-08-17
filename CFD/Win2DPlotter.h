#pragma once

#include "Win2DPlotter.g.h"

namespace winrt::CFD::implementation
{
    struct Win2DPlotter : Win2DPlotterT<Win2DPlotter>
    {
        Win2DPlotter() = default;

        int32_t MyProperty();
        void MyProperty(int32_t value);
    };
}

namespace winrt::CFD::factory_implementation
{
    struct Win2DPlotter : Win2DPlotterT<Win2DPlotter, implementation::Win2DPlotter>
    {
    };
}
