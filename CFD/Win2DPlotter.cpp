#include "pch.h"
#include "Win2DPlotter.h"
#if __has_include("Win2DPlotter.g.cpp")
#include "Win2DPlotter.g.cpp"
#endif

namespace winrt::CFD::implementation
{
    int32_t Win2DPlotter::MyProperty()
    {
        throw hresult_not_implemented();
    }

    void Win2DPlotter::MyProperty(int32_t /*value*/)
    {
        throw hresult_not_implemented();
    }
}
