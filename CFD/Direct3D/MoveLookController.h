#pragma once
#define ROTATION_GAIN 0.004f    // Sensitivity adjustment for the look controller
#define MOVEMENT_GAIN 0.1f

class MoveLookController
{
private:
    DirectX::XMFLOAT3 m_position;               
    float m_pitch, m_yaw;           

    
    bool m_moveInUse;               
    UINT32 m_movePointerID;         
    DirectX::XMFLOAT2 m_moveFirstDown;          
    DirectX::XMFLOAT2 m_movePointerPosition;   
    DirectX::XMFLOAT3 m_moveCommand;            

    
    bool m_lookInUse;               
    UINT32 m_lookPointerID;         
    DirectX::XMFLOAT2 m_lookLastPoint;          
    DirectX::XMFLOAT2 m_lookLastDelta;          

    bool m_forward, m_back;         
    bool m_left, m_right;
    bool m_up, m_down;

public: 

    MoveLookController(_In_ winrt::Windows::UI::Core::CoreWindow const& window);
    void OnPointerPressed(
        _In_ winrt::Windows::UI::Core::CoreWindow const& sender,
        _In_ winrt::Windows::UI::Core::PointerEventArgs const& args
    );
    void OnPointerMoved(
        _In_ winrt::Windows::UI::Core::CoreWindow const& sender,
        _In_ winrt::Windows::UI::Core::PointerEventArgs const& args
    );
    void OnPointerReleased(
        _In_ winrt::Windows::UI::Core::CoreWindow const& sender,
        _In_ winrt::Windows::UI::Core::PointerEventArgs const& args
    );
    void OnKeyDown(
        _In_ winrt::Windows::UI::Core::CoreWindow const& sender,
        _In_ winrt::Windows::UI::Core::KeyEventArgs const& args
    );
    void OnKeyUp(
        _In_ winrt::Windows::UI::Core::CoreWindow const& sender,
        _In_ winrt::Windows::UI::Core::KeyEventArgs const& args
    );

    //MoveLookController(_In_ winrt::Windows::UI::Core::CoreWindow const& window);

    void InitWindow(_In_ winrt::Windows::UI::Core::CoreWindow const& window);

    void Orientation(_In_ float pitch, _In_ float yaw);
    DirectX::XMFLOAT3 Orientation();

    void Position(_In_ DirectX::XMFLOAT3 pos);
    DirectX::XMFLOAT3 Position();

    void Update();

    float Pitch();
    void  Pitch(_In_ float pitch);
    float Yaw();
    void  Yaw(_In_ float yaw);

};

