// A constant buffer that stores the three basic column-major matrices for composing geometry.
cbuffer ConstantBufferChangeOnResize : register(b0)
{
    matrix projection;
};

cbuffer ConstantBufferChangesEveryFrame : register(b1)
{
    matrix view;
};

cbuffer ConstantBufferChangesEveryPrim : register(b2)
{
    matrix world;
};


struct VertexShaderInput
{
	float3 pos : POSITION;
	float3 color : COLOR0;
};

struct PixelShaderInput
{
	float4 pos : SV_POSITION;
	float3 color : COLOR0;
};


PixelShaderInput main(VertexShaderInput input)
{
	PixelShaderInput output;
	float4 pos = float4(input.pos, 1.0f);

	output.pos = mul(mul(mul(pos, world), view), projection);

	output.color = input.color;

	return output;
}