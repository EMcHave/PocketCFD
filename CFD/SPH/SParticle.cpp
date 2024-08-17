#include "pch.h"
#include "SParticle.h"

namespace SPH
{
	float2::float2(float x, float y)
	{
		this->x = x;
		this->y = y;
	}

	float float2::abs(float2 const& r)
	{
		return sqrt(r.x * r.x + r.y * r.y);
	}

	float float2::dot(float2& f2)
	{
		return x * f2.x + y * f2.y;
	}

	float2 float2::operator+(float2 f2) const
	{
		return float2(x + f2.x, y + f2.y);
	}

	float2 float2::operator-(float2 f2) const
	{
		return float2(x - f2.x, y - f2.y);
	}

	float2& float2::operator+=(float2 f2)
	{
		x += f2.x;
		y += f2.y;
		return *this;
	}

	float2& float2::operator-=(float2 f2)
	{
		x -= f2.x;
		y -= f2.y;
		return *this;
	}

	float& float2::operator[](int i)
	{
		if (i == 0) return x;
		else return y;
	}

	SParticle::SParticle(float2 r, float2 v, float R): R(R)
	{
		this->r(r);
		this->v(v);
	}

	float2 operator*(float2 f2, float f)
	{
		return float2(f * f2.x, f * f2.y);
	}

	float2 operator*(float f, float2 f2)
	{
		return float2(f * f2.x, f * f2.y);
	}

	float2 operator/(float2 f2, float f)
	{
		return float2(f2.x / f, f2.y / f);
	}

	float2 operator/(float f, float2 f2)
	{
		return float2(f2.x / f, f2.y / f);
	}

}

