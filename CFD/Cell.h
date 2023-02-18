#pragma once
#include "Node.h"
struct Cell
{
	int ID;

	double dx;
	double dy;
	double PRES_dtau;

	Node* n1;
	Node* n2;
	Node* n3;
	Node* n4;

	double p;
	double p_k;
	double xp;

	Cell(Node* n1, Node* n2,
		Node* n3, Node* n4,
		int ID, int NT,
		double dx, double dy, double PRES_dtau)
		:n1(n1), n2(n2), n3(n3), n4(n4),
		ID(ID), dx(dx), dy(dy),
		p(0), p_k(0), xp(0)
	{
		n1->vol += dx * dy / 4;
		n2->vol += dx * dy / 4;
		n3->vol += dx * dy / 4;
		n4->vol += dx * dy / 4;

		u12 = 0, u34 = 0, v14 = 0, v23 = 0;
		u14 = 0, u23 = 0, v12 = 0, v34 = 0;

		this->PRES_dtau = PRES_dtau;
	}
	~Cell() {  }

	void EvaluateIntegralsForCell(double dt, double dtau, double rho, double nu, int n);
	void XStepForCell(double dt, double dtau, double rho, double nu);
	void YStepForCell(double dt, double dtau, double rho, double nu);
	void EvaluateXp(enum IterationStep);
	void UpdatePressure(bool IsIterationsBreak);
	void UpdateAverageVelocities();
	double Divergence()
	{
		return xp / PRES_dtau;
	}
	double rot();
private:
	double u12, u34, v14, v23;
	double u14, u23, v12, v34;

};

enum IterationStep
{
	xStep,
	yStep
};