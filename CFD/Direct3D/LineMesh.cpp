#include "pch.h"
#include "LineMesh.h"
#include "ConstantBuffers.h"


LineMesh::LineMesh(bool isLine, Axis axis, winrt::com_ptr<ID3D11Device3> const& device) : MeshObject(isLine)
{
    const D3D11_INPUT_ELEMENT_DESC PNTVertexLayout[] =
    {
        { "POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0,  0, D3D11_INPUT_PER_VERTEX_DATA, 0 },
        { "COLOR",    0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 12, D3D11_INPUT_PER_VERTEX_DATA, 0 },
    };

    float3 color;
    switch (axis)
    {
    case X:
        color = float3(1.f, 0.f, 0.f);
        break;
    case Y:
        color = float3(0.f, 1.f, 0.f);
        break;
    case Z:
        color = float3(0.f, 0.f, 1.f);
        break;
    case W:
        color = float3(1.f, 1.f, 1.f);
        break;
    default:
        break;
    }

    PNTVertex cubeVertices[] =
    {
        { float3(-0.5f, 0.f, 0.f), color }, // +Y (top face)
        { float3(0.5f, 0.f, 0.f), color },
    };

    unsigned short cubeIndices[] =
    {
        0, 1
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
