#include "pch.h"
#include "UnderConstructionPage.h"
#if __has_include("UnderConstructionPage.g.cpp")
#include "UnderConstructionPage.g.cpp"
#endif

using namespace winrt;
using namespace Windows::UI::Xaml;

namespace winrt::CFD::implementation
{
    int32_t UnderConstructionPage::MyProperty()
    {
        throw hresult_not_implemented();
    }

    void UnderConstructionPage::MyProperty(int32_t /* value */)
    {
        throw hresult_not_implemented();
    }

}
