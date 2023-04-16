#pragma once
#include "pch.h"
struct Node
{
	int ID;
	bool isBoundary;
	double x;
	double y;

	double* u;
	double* v;

	double u_k;
	double xu;

	double v_k;
	double xv;

	double vol;
	double const_c;
	double integral_u;
	double integral_v;

	double aux, cux, bux, fux;
	double avx, cvx, bvx, fvx;

	double auy, cuy, buy, fuy;
	double avy, cvy, bvy, fvy;

	Node(double x, double y, int ID, int NT, double dt, double dtau)
	{
		memset(this, 0, sizeof(Node));
		this->const_c = 1 / (2 * dt) + 1 / dtau;
		this->x = x;
		this->y = y;
		this->ID = ID;
		this->isBoundary = false;

		u = new double[NT];
		v = new double[NT];
		for (int i = 0; i < NT; i++)
		{
			u[i] = 0; 
			v[i] = 0;
		}
	}
	~Node()
	{
		delete[] u;
		delete[] v;
	}

	void ResetNode(int n)
	{
		integral_u = 0;
		integral_v = 0;
		aux = 0, bux = 0, cux = 0; fux = 0;
		avx = 0, bvx = 0, cvx = 0; fvx = 0;
		auy = 0, buy = 0, cuy = 0; fuy = 0;
		avy = 0, bvy = 0, cvy = 0; fvy = 0;

		//u[n] = u[n - 1];
		//v[n] = v[n - 1];

		cux += const_c;
		cvx += const_c;
		cuy += const_c;
		cvy += const_c;
	}
	void UpdateVelocities(bool IsIterationsBreak, int n)
	{
		if (IsIterationsBreak)
		{
			u[n] = u_k;
			v[n] = v_k;
		}
		else
		{
			u_k += xu;
			v_k += xv;
		}
	}
	void EvaluateXiuv(double dtau, double dt, int n)
	{
		xu = dtau * (-(u_k - u[n - 1]) / dt + 1 / vol * integral_u) * !isBoundary;
		xv = dtau * (-(v_k - v[n - 1]) / dt + 1 / vol * integral_v) * !isBoundary;
		fux += xu / dtau;
		fvx += xv / dtau;
	}

	double AbsVelocity(int n)
	{
		return sqrt(u[n] * u[n] + v[n] * v[n]);
	}
};

