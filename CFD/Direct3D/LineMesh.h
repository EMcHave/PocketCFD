#pragma once
#include "MeshObject.h"
#include "Line.h"
class LineMesh : public MeshObject
{
public:
	LineMesh(_In_opt_ bool isLine,Axis axis, _In_ winrt::com_ptr<ID3D11Device3> const& device);
};

