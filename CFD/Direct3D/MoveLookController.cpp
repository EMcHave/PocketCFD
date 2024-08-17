#include "pch.h"
#include "MoveLookController.h"

MoveLookController::MoveLookController(winrt::Windows::UI::Core::CoreWindow const& window)
{
    InitWindow(window);

}

void MoveLookController::OnPointerPressed(winrt::Windows::UI::Core::CoreWindow const& sender, winrt::Windows::UI::Core::PointerEventArgs const& args)
{
    UINT32 pointerID = args.CurrentPoint().PointerId();
    DirectX::XMFLOAT2 position = DirectX::XMFLOAT2(args.CurrentPoint().Position().X, args.CurrentPoint().Position().Y);

    auto device = args.CurrentPoint().PointerDevice();
    auto deviceType = device.PointerDeviceType();
    if (deviceType == winrt::Windows::Devices::Input::PointerDeviceType::Mouse)
    {
        // Action, Jump, or Fire
    }

    // Check  if this pointer is in the move control.
    // Change the values  to percentages of the preferred screen resolution.
    // You can set the x value to <preferred resolution> * <percentage of width>
    // for example, ( position.x < (screenResolution.x * 0.15) ).

    if ((position.x < 300 && position.y > 380) && (deviceType != winrt::Windows::Devices::Input::PointerDeviceType::Mouse))
    {
        if (!m_moveInUse) 
        {
            m_moveFirstDown = position;                 
            m_movePointerPosition = position;
            m_movePointerID = pointerID;                
            m_moveInUse = TRUE;
        }
    }
    else 
    {
        if (!m_lookInUse) 
        {
            m_lookLastPoint = position;                         
            m_lookPointerID = args.CurrentPoint().PointerId();  
            m_lookLastDelta.x = m_lookLastDelta.y = 0;          
            m_lookInUse = TRUE;
        }
    }
}

void MoveLookController::OnPointerMoved(winrt::Windows::UI::Core::CoreWindow const& sender, winrt::Windows::UI::Core::PointerEventArgs const& args)
{
    UINT32 pointerID = args.CurrentPoint().PointerId();
    DirectX::XMFLOAT2 position = DirectX::XMFLOAT2(args.CurrentPoint().Position().X, args.CurrentPoint().Position().Y);

    // Decide which control this pointer is operating.
    if (pointerID == m_movePointerID)           // This is the move pointer.
    {
        // Move control
        m_movePointerPosition = position;       // Save the current position.

    }
    else if (pointerID == m_lookPointerID)      // This is the look pointer.
    {
        // Look control

        DirectX::XMFLOAT2 pointerDelta;
        pointerDelta.x = position.x - m_lookLastPoint.x;        // How far did pointer move
        pointerDelta.y = position.y - m_lookLastPoint.y;

        DirectX::XMFLOAT2 rotationDelta;
        rotationDelta.x = pointerDelta.x * ROTATION_GAIN;   // Scale for control sensitivity.
        rotationDelta.y = pointerDelta.y * ROTATION_GAIN;

        m_lookLastPoint = position;                     // Save for the next time through.

        // Update our orientation based on the command.
        m_pitch -= rotationDelta.y;                     // Mouse y increases down, but pitch increases up.
        m_yaw -= rotationDelta.x;                       // Yaw is defined as CCW around the y-axis.

        // Limit the pitch to straight up or straight down.
        m_pitch = (float)__max(-DirectX::XM_PI / 2.0f, m_pitch);
        m_pitch = (float)__min(+DirectX::XM_PI / 2.0f, m_pitch);
    }
}

void MoveLookController::OnPointerReleased(winrt::Windows::UI::Core::CoreWindow const& sender, winrt::Windows::UI::Core::PointerEventArgs const& args)
{
    UINT32 pointerID = args.CurrentPoint().PointerId();
    DirectX::XMFLOAT2 position = DirectX::XMFLOAT2(args.CurrentPoint().Position().X, args.CurrentPoint().Position().Y);


    if (pointerID == m_movePointerID)    // This was the move pointer.
    {
        m_moveInUse = FALSE;
        m_movePointerID = 0;
    }
    else if (pointerID == m_lookPointerID) // This was the look pointer.
    {
        m_lookInUse = FALSE;
        m_lookPointerID = 0;
    }
}

void MoveLookController::OnKeyDown(winrt::Windows::UI::Core::CoreWindow const& sender, winrt::Windows::UI::Core::KeyEventArgs const& args)
{
    using namespace winrt::Windows::System;
    winrt::Windows::System::VirtualKey Key;
    Key = args.VirtualKey();

    // Figure out the command from the keyboard.
    if (Key == VirtualKey::W)     // Forward
        m_forward = true;
    if (Key == VirtualKey::S)     // Back
        m_back = true;
    if (Key == VirtualKey::A)     // Left
        m_left = true;
    if (Key == VirtualKey::D)     // Right
        m_right = true;
    if (Key == VirtualKey::X)
        m_up = true;
    if (Key == VirtualKey::Space)
        m_down = true;
}

void MoveLookController::OnKeyUp(winrt::Windows::UI::Core::CoreWindow const& sender, winrt::Windows::UI::Core::KeyEventArgs const& args)
{
    using namespace winrt::Windows::System;
    winrt::Windows::System::VirtualKey Key;
    Key = args.VirtualKey();

    // Figure out the command from the keyboard.
    if (Key == VirtualKey::W)     // forward
        m_forward = false;
    if (Key == VirtualKey::S)     // back
        m_back = false;
    if (Key == VirtualKey::A)     // left
        m_left = false;
    if (Key == VirtualKey::D)     // right
        m_right = false;
    if (Key == VirtualKey::X)
        m_up = false;
    if (Key == VirtualKey::Space)
        m_down = false;
}

void MoveLookController::InitWindow(winrt::Windows::UI::Core::CoreWindow const& window)
{
    m_moveInUse = false;
    m_lookInUse = false;
    m_movePointerID = 0;
    m_lookPointerID = 0;
    m_moveCommand = DirectX::XMFLOAT3(0.0f, 0.0f, 0.0f);
    m_forward = false;
    m_back = false;
    m_left = false;
    m_right = false;
    m_up = false;
    m_down = false;

    window.PointerPressed({ this, &MoveLookController::OnPointerPressed });

    window.PointerMoved({ this, &MoveLookController::OnPointerMoved });

    window.PointerReleased({ this, &MoveLookController::OnPointerReleased });

    window.KeyDown({ this, &MoveLookController::OnKeyDown });

    window.KeyUp({ this, &MoveLookController::OnKeyUp });

    m_moveInUse = FALSE;                // No pointer is in the Move control.
    m_movePointerID = 0;

    m_lookInUse = FALSE;                // No pointer is in the Look control.
    m_lookPointerID = 0;

    //  Need to init this as it is reset every frame.
    m_moveCommand = DirectX::XMFLOAT3(0.0f, 0.0f, 0.0f);

    m_position = DirectX::XMFLOAT3(-1.6f, 7.9f, 9.6f);
    Orientation(-0.65, -3.61);
}

void MoveLookController::Orientation(float pitch, float yaw)
{
    m_pitch = pitch;
    m_yaw = yaw;
}

DirectX::XMFLOAT3 MoveLookController::Orientation()
{
    float y = sinf(m_pitch);        // Vertical
    float r = cosf(m_pitch);        // In the plane
    float z = r * cosf(m_yaw);        // Fwd-back
    float x = r * sinf(m_yaw);        // Left-right
    DirectX::XMFLOAT3 result(x, y, z);
    result.x += m_position.x;
    result.y += m_position.y;
    result.z += m_position.z;


    return result;
}

void MoveLookController::Position(DirectX::XMFLOAT3 pos)
{
    m_position = pos;
}

DirectX::XMFLOAT3 MoveLookController::Position()
{
    return m_position;
}

void MoveLookController::Update()
{
    if (m_moveInUse)
    {
        DirectX::XMFLOAT2 pointerDelta(m_movePointerPosition);
        pointerDelta.x -= m_moveFirstDown.x;
        pointerDelta.y -= m_moveFirstDown.y;

        // Figure out the command from the touch-based virtual joystick.
        if (pointerDelta.x > 16.0f)      // Leave 32 pixel-wide dead spot for being still.
            m_moveCommand.x = 1.0f;
        else
            if (pointerDelta.x < -16.0f)
                m_moveCommand.x = -1.0f;

        if (pointerDelta.y > 16.0f)      // Joystick y is up, so change sign.
            m_moveCommand.y = -1.0f;
        else
            if (pointerDelta.y < -16.0f)
                m_moveCommand.y = 1.0f;
    }

    // Poll our state bits that are set by the keyboard input events.
    if (m_forward)
        m_moveCommand.y += 1.0f;
    if (m_back)
        m_moveCommand.y -= 1.0f;

    if (m_left)
        m_moveCommand.x -= 1.0f;
    if (m_right)
        m_moveCommand.x += 1.0f;

    if (m_up)
        m_moveCommand.z += 1.0f;
    if (m_down)
        m_moveCommand.z -= 1.0f;

    // Make sure that 45 degree cases are not faster.
    DirectX::XMFLOAT3 command = m_moveCommand;
    DirectX::XMVECTOR vector;
    vector = DirectX::XMLoadFloat3(&command);

    if (fabsf(command.x) > 0.1f || fabsf(command.y) > 0.1f || fabsf(command.z) > 0.1f)
    {
        vector = DirectX::XMVector3Normalize(vector);
        DirectX::XMStoreFloat3(&command, vector);
    }


    // Rotate command to align with our direction (world coordinates).
    DirectX::XMFLOAT3 wCommand;
    wCommand.x = command.x * cosf(m_yaw) - command.y * sinf(m_yaw);
    wCommand.y = command.x * sinf(m_yaw) + command.y * cosf(m_yaw);
    wCommand.z = command.z;

    // Scale for sensitivity adjustment.
    wCommand.x = wCommand.x * MOVEMENT_GAIN;
    wCommand.y = wCommand.y * MOVEMENT_GAIN;
    wCommand.z = wCommand.z * MOVEMENT_GAIN;

    // Our velocity is based on the command.
    // Also note that y is the up-down axis. 
    DirectX::XMFLOAT3 Velocity;
    Velocity.x = -wCommand.x;
    Velocity.z = wCommand.y;
    Velocity.y = wCommand.z;


    // Integrate
    m_position.x += Velocity.x;
    m_position.y += Velocity.y;
    m_position.z += Velocity.z;

    // Clear movement input accumulator for use during the next frame.
    m_moveCommand = DirectX::XMFLOAT3(0.0f, 0.0f, 0.0f);
}

float MoveLookController::Pitch()
{
    return m_pitch;
}

void MoveLookController::Pitch(_In_ float pitch)
{
    m_pitch = pitch;
}

float MoveLookController::Yaw()
{
    return m_yaw;
}

void MoveLookController::Yaw(_In_ float yaw)
{
    m_yaw = yaw;
}