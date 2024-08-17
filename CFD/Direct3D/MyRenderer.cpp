#include "pch.h"
#include "MyRenderer.h"
#include "Cube.h"
#include "CubeMesh.h"
#include "LineMesh.h"


using namespace DirectX;

MyRenderer::MyRenderer(std::shared_ptr<DX::DeviceResources> const& deviceResources) :
    m_deviceResources(deviceResources)
{
    //CreateDeviceDependentResources();
    //CreateWindowSizeDependentResources();
}

void MyRenderer::CreateDeviceDependentResources(_In_ std::shared_ptr<LogicClass> logic)
{
    m_objects = logic->RenderObjects();

    m_lines = std::vector<std::shared_ptr<Line>>();
    m_lines.push_back(std::make_shared<Line>(Axis::X));
    m_lines.push_back(std::make_shared<Line>(Axis::Y));
    m_lines.push_back(std::make_shared<Line>(Axis::Z));



    m_game = logic;

    using namespace winrt::Windows::Foundation::Numerics;
    
    auto d3dDevice = m_deviceResources->GetD3DDevice();

    const D3D11_INPUT_ELEMENT_DESC PNTVertexLayout[] =
    {
        { "POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0,  0, D3D11_INPUT_PER_VERTEX_DATA, 0 },
        { "COLOR",    0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 12, D3D11_INPUT_PER_VERTEX_DATA, 0 },
    };

    D3D11_BUFFER_DESC bd;
    ZeroMemory(&bd, sizeof(bd));

    bd.Usage = D3D11_USAGE_DEFAULT;
    bd.BindFlags = D3D11_BIND_CONSTANT_BUFFER;
    bd.CPUAccessFlags = 0;
    bd.ByteWidth = (sizeof(ConstantBufferChangeOnResize) + 15) / 16 * 16;
    m_constantBufferChangeOnResize = nullptr;
    winrt::check_hresult(
        d3dDevice->CreateBuffer(&bd, nullptr, m_constantBufferChangeOnResize.put())
    );

    bd.ByteWidth = (sizeof(ConstantBufferChangesEveryFrame) + 15) / 16 * 16;
    m_constantBufferChangesEveryFrame = nullptr;
    winrt::check_hresult(
        d3dDevice->CreateBuffer(&bd, nullptr, m_constantBufferChangesEveryFrame.put())
    );

    bd.ByteWidth = (sizeof(ConstantBufferChangesEveryPrim) + 15) / 16 * 16;
    m_constantBufferChangesEveryPrim = nullptr;
    winrt::check_hresult(
        d3dDevice->CreateBuffer(&bd, nullptr, m_constantBufferChangesEveryPrim.put())
    );

    BasicLoader loader{ d3dDevice };

    uint32_t numElements = ARRAYSIZE(PNTVertexLayout);
    loader.LoadShader(L"VertexShader.cso", PNTVertexLayout, numElements, m_vertexShader.put(), m_inputLayout.put());
    loader.LoadShader(L"PixelShader.cso", m_pixelShader.put());

    auto cubeMesh = std::make_shared<CubeMesh>(false, d3dDevice);



    
    for (auto&& object : m_objects)
        object->Mesh(cubeMesh);

    m_lines.at(0)->Mesh(std::make_shared<LineMesh>(true, Axis::X, d3dDevice));
    m_lines.at(1)->Mesh(std::make_shared<LineMesh>(true, Axis::Y, d3dDevice));
    m_lines.at(2)->Mesh(std::make_shared<LineMesh>(true, Axis::Z, d3dDevice));
    for (int i = 3; i < m_lines.size(); i++)
        m_lines.at(i)->Mesh(std::make_shared<LineMesh>(true, Axis::W, d3dDevice));
    
}

void MyRenderer::CreateWindowSizeDependentResources()
{
    degree = 0;
    auto d3dContext = m_deviceResources->GetD3DDeviceContext();
    auto outputSize = m_deviceResources->GetOutputSize();
    float aspectRatio = outputSize.Width / outputSize.Height;
    float fovAngleY = 70.0f * XM_PI / 180.0f;

    if (aspectRatio < 1.0f)
    {
        fovAngleY *= 2.0f;
    }

    XMMATRIX perspectiveMatrix = XMMatrixPerspectiveFovLH(
        fovAngleY,
        aspectRatio,
        0.5f,
        10.0f
    );

    m_game->GameCamera().SetProjParams(
        fovAngleY, aspectRatio,
        0.01f,
        100.0f
    );

    XMFLOAT4X4 orientation = m_deviceResources->GetOrientationTransform3D();

    XMMATRIX orientationMatrix = XMLoadFloat4x4(&orientation);

    ConstantBufferChangeOnResize changesOnResize;
    XMStoreFloat4x4(&changesOnResize.projection,
        XMMatrixMultiply(
            XMMatrixTranspose(m_game->GameCamera().Projection()) /*XMMatrixTranspose(perspectiveMatrix * orientationMatrix)*/,
            XMMatrixTranspose(orientationMatrix))
            );


    static const XMVECTORF32 eye = { 0.0f, 2.0f, 1.f, 0.0f };
    static const XMVECTORF32 at = { 0.0f, 0.0f, 0.0f, 0.0f };
    static const XMVECTORF32 up = { 0.0f, 1.0f, 0.0f, 0.0f };

    view = XMMatrixTranspose(XMMatrixLookAtLH(eye, at, up));

    d3dContext->UpdateSubresource(
        m_constantBufferChangeOnResize.get(),
        0,
        nullptr,
        &changesOnResize,
        0,
        0
    );
}

void MyRenderer::ReleaseDeviceDependentResources()
{
}


void MyRenderer::Render()
{
    degree += 0.1;

    auto m_d3dDeviceContext = m_deviceResources->GetD3DDeviceContext();
    const float clearColor[4] = { 0.071f, 0.04f, 0.561f, 1.0f };
    //const float clearColor[4] = { 1.f, 1.f, 1.f, 0.0f };


    ID3D11RenderTargetView* const targets[1] = { m_deviceResources->GetBackBufferRenderTargetView() };

    m_d3dDeviceContext->OMSetRenderTargets(
        1,
        targets,
        m_deviceResources->GetDepthStencilView()
    );
    m_d3dDeviceContext->ClearDepthStencilView(
        m_deviceResources->GetDepthStencilView(),
        D3D11_CLEAR_DEPTH,
        1.0f,
        0
    );
    m_d3dDeviceContext->ClearRenderTargetView(
        m_deviceResources->GetBackBufferRenderTargetView(),
        clearColor
    );


    ConstantBufferChangesEveryFrame constantBufferChangesEveryFrameValue;
    XMStoreFloat4x4(
        &constantBufferChangesEveryFrameValue.view,
        XMMatrixTranspose(m_game->GameCamera().View())
    );
    m_d3dDeviceContext->UpdateSubresource(
        m_constantBufferChangesEveryFrame.get(),
        0,
        nullptr,
        &constantBufferChangesEveryFrameValue,
        0,
        0
    );



    m_d3dDeviceContext->IASetInputLayout(m_inputLayout.get());

    // Set the vertex and pixel shader stage state.
    ID3D11Buffer* constantBufferChangeOnResize{ m_constantBufferChangeOnResize.get() };
    m_d3dDeviceContext->VSSetConstantBuffers(0, 1, &constantBufferChangeOnResize);
    ID3D11Buffer* constantBufferChangesEveryFrame{ m_constantBufferChangesEveryFrame.get() };
    m_d3dDeviceContext->VSSetConstantBuffers(1, 1, &constantBufferChangesEveryFrame);
    ID3D11Buffer* constantBufferChangesEveryPrim{ m_constantBufferChangesEveryPrim.get() };
    m_d3dDeviceContext->VSSetConstantBuffers(2, 1, &constantBufferChangesEveryPrim);

    m_d3dDeviceContext->PSSetConstantBuffers(1, 1, &constantBufferChangesEveryFrame);
    m_d3dDeviceContext->PSSetConstantBuffers(2, 1, &constantBufferChangesEveryPrim);

    m_d3dDeviceContext->VSSetShader(
        m_vertexShader.get(),
        nullptr,
        0
    );

    m_d3dDeviceContext->PSSetShader(
        m_pixelShader.get(),
        nullptr,
        0
    );
    
    for (auto&& object : m_objects)
        object->Render(m_d3dDeviceContext, m_constantBufferChangesEveryPrim.get());
    for (auto&& object : m_lines)
        object->Render(m_d3dDeviceContext, m_constantBufferChangesEveryPrim.get());
    //for (int i = 3; i < m_lines.size(); i++)
        //m_lines[i]->Render(m_d3dDeviceContext, m_constantBufferChangesEveryPrim.get());
}
