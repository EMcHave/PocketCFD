#include "pch.h"
#include "MainPage.h"
#include "MainPage.g.cpp"
#include <winrt/Windows.UI.Xaml.Interop.h>

using namespace winrt;
using namespace Windows::UI::Xaml;
using namespace Microsoft::UI::Xaml::Controls;

namespace winrt::CFD::implementation
{
    MainPage::MainPage()
    {
        ApplicationViewTitleBar appTitleBar =
            ApplicationView::GetForCurrentView().TitleBar();
        //appTitleBar.ButtonBackgroundColor() = Windows::UI::Colors::Transparent();

        //auto coreTitleBar = winrt::Windows::ApplicationModel::Core::CoreApplication::GetCurrentView().TitleBar();
        //coreTitleBar.ExtendViewIntoTitleBar(true);
        //Window::Current().SetTitleBar(AppTitleBar());
    }
}


void winrt::CFD::implementation::MainPage::Page_Loaded(IInspectable const& sender, RoutedEventArgs const& e)
{

}

void winrt::CFD::implementation::MainPage::nvView_ItemInvoked(NavigationView const& sender, NavigationViewItemInvokedEventArgs const& args)
{
    Interop::TypeName pageTypeName;
    //Interop::TypeName pageTypeName2 = xaml_typename<CFD::CFDPage>();
    pageTypeName.Name = L"CFD." + unbox_value<hstring>(args.InvokedItemContainer().Tag());
    pageTypeName.Kind = Interop::TypeKind::Primitive;
    contentFrame().Navigate(pageTypeName, nullptr);
}


void winrt::CFD::implementation::MainPage::nvView_Loaded(IInspectable const& sender, RoutedEventArgs const& e)
{
    nvView().SelectedItem(nvView().MenuItems().GetAt(0));
    m_pages.push_back(std::make_pair<std::wstring, Windows::UI::Xaml::Interop::TypeName>
        (L"CFDPage", winrt::xaml_typename<CFD::CFDPage>()));
    m_pages.push_back(std::make_pair<std::wstring, Windows::UI::Xaml::Interop::TypeName>
        (L"UnderConstructionPage", winrt::xaml_typename<CFD::UnderConstructionPage>()));

    nvView_Navigate(L"CFDPage", Windows::UI::Xaml::Media::Animation::EntranceNavigationTransitionInfo());
}


void winrt::CFD::implementation::MainPage::nvView_SelectionChanged(NavigationView const& sender, NavigationViewSelectionChangedEventArgs const& args)
{

}

void winrt::CFD::implementation::MainPage::nvView_Navigate(
    std::wstring navItemTag,
    Windows::UI::Xaml::Media::Animation::NavigationTransitionInfo const& transitionInfo)
{
    Windows::UI::Xaml::Interop::TypeName pageTypeName;
    for (auto&& eachPage : m_pages)
    {
        if (eachPage.first == navItemTag)
        {
            pageTypeName = eachPage.second;
            break;
        }
    }

    Windows::UI::Xaml::Interop::TypeName preNavPageType =
        contentFrame().CurrentSourcePageType();

    
    if (pageTypeName.Name != L"" && preNavPageType.Name != pageTypeName.Name)
    {
        contentFrame().Navigate(pageTypeName, nullptr, transitionInfo);
        if (navItemTag == L"CFDPage")
            nvView().Header(box_value(L"CFD"));
    }
}