#pragma once

#include "UnderConstructionPage.g.h"

namespace winrt::CFD::implementation
{
    struct UnderConstructionPage : UnderConstructionPageT<UnderConstructionPage>
    {
        UnderConstructionPage() 
        {
            // Xaml objects should not call InitializeComponent during construction.
            // See https://github.com/microsoft/cppwinrt/tree/master/nuget#initializecomponent
        }

        int32_t MyProperty();
        void MyProperty(int32_t value);

    
    };
}

namespace winrt::CFD::factory_implementation
{
    struct UnderConstructionPage : UnderConstructionPageT<UnderConstructionPage, implementation::UnderConstructionPage>
    {
    };
}
