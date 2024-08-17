#include "pch.h"
#include "CubeMesh.h"
#include "ConstantBuffers.h"

using namespace winrt::Windows::Foundation::Numerics;

CubeMesh::CubeMesh(_In_opt_ bool isLine, _In_ winrt::com_ptr<ID3D11Device3> const& device) : MeshObject(_In_opt_ isLine)
{
    using namespace DirectX;
    const D3D11_INPUT_ELEMENT_DESC PNTVertexLayout[] =
    {
        { "POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0,  0, D3D11_INPUT_PER_VERTEX_DATA, 0 },
        { "COLOR",    0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 12, D3D11_INPUT_PER_VERTEX_DATA, 0 },
    };

    PNTVertex cubeVertices[] =
    {
        { float3(-0.5f, 0.5f, -0.5f), float3(0.0f, 1.0f, 0.0f) }, // +Y (top face)
        { float3(0.5f, 0.5f, -0.5f), float3(1.0f, 1.0f, 0.0f) },
        { float3(0.5f, 0.5f,  0.5f), float3(1.0f, 1.0f, 1.0f) },
        { float3(-0.5f, 0.5f,  0.5f), float3(0.0f, 1.0f, 1.0f) },

        { float3(-0.5f, -0.5f,  0.5f), float3(0.0f, 0.0f, 1.0f) }, // -Y (bottom face)
        { float3(0.5f, -0.5f,  0.5f), float3(1.0f, 0.0f, 1.0f) },
        { float3(0.5f, -0.5f, -0.5f), float3(1.0f, 0.0f, 0.0f) },
        { float3(-0.5f, -0.5f, -0.5f), float3(0.0f, 0.0f, 0.0f) },
    };

    unsigned short cubeIndices[] =
    {
        0, 1, 2,
        0, 2, 3,

        4, 5, 6,
        4, 6, 7,

        3, 2, 5,
        3, 5, 4,

        2, 1, 6,
        2, 6, 5,

        1, 7, 6,
        1, 0, 7,

        0, 3, 4,
        0, 4, 7
    };

    m_vertexCount = ARRAYSIZE(cubeVertices);
    m_indexCount = ARRAYSIZE(cubeIndices);

    D3D11_BUFFER_DESC vertexBufferDesc = { 0 };
    D3D11_SUBRESOURCE_DATA vertexBufferData = { 0 };

    //ZeroMemory(&vertexBufferDesc, sizeof(vertexBufferDesc));
    vertexBufferDesc.ByteWidth = sizeof(PNTVertex) * m_vertexCount;
    vertexBufferDesc.Usage = D3D11_USAGE_DEFAULT;
    vertexBufferDesc.BindFlags = D3D11_BIND_VERTEX_BUFFER;
    vertexBufferDesc.CPUAccessFlags = 0;

    vertexBufferData.pSysMem = cubeVertices;
    vertexBufferData.SysMemPitch = 0;
    vertexBufferData.SysMemSlicePitch = 0;


    winrt::check_hresult(
        device->CreateBuffer(
            &vertexBufferDesc,
            &vertexBufferData,
            m_vertexBuffer.put()
        )
    );

    D3D11_BUFFER_DESC indexBufferDesc = { 0 };
    //ZeroMemory(&indexBufferDesc, sizeof(indexBufferDesc));
    indexBufferDesc.ByteWidth = sizeof(unsigned short) * m_indexCount;
    indexBufferDesc.Usage = D3D11_USAGE_DEFAULT;
    indexBufferDesc.BindFlags = D3D11_BIND_INDEX_BUFFER;
    indexBufferDesc.CPUAccessFlags = 0;

    D3D11_SUBRESOURCE_DATA indexBufferData = { 0 };
    indexBufferData.pSysMem = cubeIndices;
    indexBufferData.SysMemPitch = 0;
    indexBufferData.SysMemSlicePitch = 0;


    winrt::check_hresult(
        device->CreateBuffer(
            &indexBufferDesc,
            &indexBufferData,
            m_indexBuffer.put()
        )
    );
}