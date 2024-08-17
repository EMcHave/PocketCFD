#include "pch.h"
#include "GeometricObject.h"

using namespace DirectX;

GeometricObject::GeometricObject()
{
    m_position = XMFLOAT3(0.0f, 0.0f, 0.0f);
    m_velocity = XMFLOAT3(0.0f, 0.0f, 0.0f);
    m_angularVelocity = XMFLOAT3(0.0f, 0.0f, 0.0f);
    m_defaultXAxis = XMFLOAT3(1.0f, 0.0f, 0.0f);
    m_defaultYAxis = XMFLOAT3(0.0f, 1.0f, 0.0f);
    m_defaultZAxis = XMFLOAT3(0.0f, 0.0f, 1.0f);
    m_quaternion = XMFLOAT4(0.f, 0.f, 0.f, 1.f);
    XMStoreFloat4x4(&m_modelMatrix, XMMatrixIdentity());
}

void GeometricObject::Render(
    _In_ ID3D11DeviceContext* context,
    _In_ ID3D11Buffer* primitiveConstantBuffer
)
{

    ConstantBufferChangesEveryPrim constantBuffer;

    XMStoreFloat4x4(
        &constantBuffer.world,
        XMMatrixTranspose(ModelMatrix())
    );

    context->UpdateSubresource(primitiveConstantBuffer, 0, nullptr, &constantBuffer, 0, 0);

    m_mesh->Render(context);
}
