#pragma once

#include "pch.h"

using namespace winrt::Windows::Foundation::Numerics;

struct PNTVertex
{
    float3 pos;   // Position
    float3 color; // Color
};

struct ConstantBufferChangeOnResize
{
    DirectX::XMFLOAT4X4 projection;
};

struct ConstantBufferChangesEveryFrame
{
    DirectX::XMFLOAT4X4 view;
};

struct ConstantBufferChangesEveryPrim
{
    DirectX::XMFLOAT4X4 world;
};