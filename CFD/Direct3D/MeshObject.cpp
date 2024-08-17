#include "pch.h"
#include "MeshObject.h"

#include "pch.h"
#include "MeshObject.h"
#include "ConstantBuffers.h"

MeshObject::MeshObject(_In_opt_ bool isLine) :
    m_vertexCount(0),
    m_indexCount(0),
    m_isLine(isLine)
{
}

void MeshObject::Render(_In_ ID3D11DeviceContext* context)
{
    uint32_t stride{ sizeof(PNTVertex) };
    uint32_t offset{ 0 };

    ID3D11Buffer* vertexBuffer{ m_vertexBuffer.get() };
    context->IASetVertexBuffers(0, 1, &vertexBuffer, &stride, &offset);
    context->IASetIndexBuffer(m_indexBuffer.get(), DXGI_FORMAT_R16_UINT, 0);
    if (m_isLine)
        context->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_LINELIST);
    else
        context->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
    context->DrawIndexed(m_indexCount, 0, 0);
}