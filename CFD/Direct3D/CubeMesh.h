#pragma once

#include "MeshObject.h"

class CubeMesh : public MeshObject
{
public:
	CubeMesh(_In_opt_ bool isLine, _In_ winrt::com_ptr<ID3D11Device3> const& device);
};

