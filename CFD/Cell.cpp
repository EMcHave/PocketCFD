#include "pch.h"
#include "Cell.h"

void Cell::EvaluateIntegralsForCell(double dt, double dtau, double rho, double nu, int n)
{
	double Fux12, Fux34, Fuy14, Fuy23;
	double Fvx12, Fvx34, Fvy14, Fvy23;

	u12 = (n1->u_k + n2->u_k) / 2;
	u34 = (n3->u_k + n4->u_k) / 2;
	v14 = (n1->v_k + n4->v_k) / 2;
	v23 = (n2->v_k + n3->v_k) / 2;

	u14 = (n1->u_k + n4->u_k) / 2;
	u23 = (n2->u_k + n3->u_k) / 2;
	v12 = (n1->v_k + n2->v_k) / 2;
	v34 = (n3->v_k + n4->v_k) / 2;

	Fux12 = dy / 2 * (-p_k / rho + nu * (n2->u_k - n1->u_k) / dx
		- (abs(u12) + u12) / 2 * n1->u_k - (u12 - abs(u12)) / 2 * n2->u_k);
	Fux34 = dy / 2 * (-p_k / rho + nu * (n3->u_k - n4->u_k) / dx
		- (abs(u34) + u34) / 2 * n4->u_k - (u34 - abs(u34)) / 2 * n3->u_k);
	Fuy14 = dx / 2 * (nu * (n4->u_k - n1->u_k) / dy
		- (abs(v14) + v14) / 2 * n1->u_k - (v14 - abs(v14)) / 2 * n4->u_k);
	Fuy23 = dx / 2 * (nu * (n3->u_k - n2->u_k) / dy
		- (abs(v23) + v23) / 2 * n2->u_k - (v23 - abs(v23)) / 2 * n3->u_k);

	Fvx12 = dy / 2 * (nu * (n2->v_k - n1->v_k) / dx
		- (abs(u12) + u12) / 2 * n1->v_k - (u12 - abs(u12)) / 2 * n2->v_k);
	Fvx34 = dy / 2 * (nu * (n3->v_k - n4->v_k) / dx
		- (abs(u34) + u34) / 2 * n4->v_k - (u34 - abs(u34)) / 2 * n3->v_k);
	Fvy14 = dx / 2 * (-p_k / rho + nu * (n4->v_k - n1->v_k) / dy
		- (abs(v14) + v14) / 2 * n1->v_k - (v14 - abs(v14)) / 2 * n4->v_k);
	Fvy23 = dx / 2 * (-p_k / rho + nu * (n3->v_k - n2->v_k) / dy
		- (abs(v23) + v23) / 2 * n2->v_k - (v23 - abs(v23)) / 2 * n3->v_k);

	n1->integral_u += Fux12 + Fuy14;
	n2->integral_u += -Fux12 + Fuy23;
	n3->integral_u += -Fux34 - Fuy23;
	n4->integral_u += Fux34 - Fuy14;

	n1->integral_v += Fvx12 + Fvy14;
	n2->integral_v += -Fvx12 + Fvy23;
	n3->integral_v += -Fvx34 - Fvy23;
	n4->integral_v += Fvx34 - Fvy14;

	//this->xp = -PRES_dtau * ((u23 - u14) / dx + (v34 - v12) / dy) * !isFreeOutlet;
	this->xp = -PRES_dtau * ((u23 - u14) / dx + (v34 - v12) / dy) * !isFreeOutlet;
}

void Cell::XStepForCell(double dt, double dtau, double rho, double nu)
{
	double mult1 = 1 / n1->vol * dy / 2;
	double mult2 = 1 / n2->vol * dy / 2;
	double mult3 = 1 / n3->vol * dy / 2;
	double mult4 = 1 / n4->vol * dy / 2;

	double c12 = (((abs(u12) + u12)) / 2 + nu * 1 / dx);
	double b12 = (-((u12 - abs(u12))) / 2 + nu * 1 / dx);

	double c34 = (((abs(u34) + u34)) / 2 + nu * 1 / dx);
	double b34 = (-((u34 - abs(u34))) / 2 + nu * 1 / dx);

	n1->cux += (c12 + PRES_dtau / rho / dx) * mult1;
	n1->bux += (b12 + PRES_dtau / rho / dx) * mult1 * !n1->isBoundary;
	n1->fux += -mult1 * xp / rho * !n1->isBoundary;

	n4->cux += (c34 + PRES_dtau / rho / dx) * mult4;
	n4->bux += (b34 + PRES_dtau / rho / dx) * mult4 * !n4->isBoundary;
	n4->fux += -mult4 * xp / rho * !n4->isBoundary;

	n2->aux += (c12 + PRES_dtau / rho / dx) * mult2 * !n2->isBoundary;
	n2->cux += (b12 + PRES_dtau / rho / dx) * mult2;
	n2->fux += mult2 * xp / rho * !n2->isBoundary;

	n3->aux += (c34 + PRES_dtau / rho / dx) * mult3 * !n3->isBoundary;
	n3->cux += (b34 + PRES_dtau / rho / dx) * mult3;
	n3->fux += mult3 * xp / rho * !n3->isBoundary;


	n1->cvx += mult1 * c12;
	n1->bvx += mult2 * b12 * !n1->isBoundary;

	n4->cvx += mult4 * c34;
	n4->bvx += mult4 * b34 * !n4->isBoundary;

	n2->avx += mult2 * c12 * !n2->isBoundary;
	n2->cvx += mult2 * b12;

	n3->avx += mult3 * c34 * !n3->isBoundary;
	n3->cvx += mult3 * b34;
}

void Cell::YStepForCell(double dt, double dtau, double rho, double nu)
{
	double mult1 = 1 / n1->vol * dx / 2;
	double mult2 = 1 / n2->vol * dx / 2;
	double mult3 = 1 / n3->vol * dx / 2;
	double mult4 = 1 / n4->vol * dx / 2;

	double c14 = (((abs(v14) + v14)) / 2 + nu * 1 / dy);
	double b14 = (-((v14 - abs(v14))) / 2 + nu * 1 / dy);
	double c23 = (((abs(v23) + v23)) / 2 + nu * 1 / dy);
	double b23 = (-((v23 - abs(v23))) / 2 + nu * 1 / dy);

	n1->cuy += mult1 * c14;
	n1->buy += mult1 * b14 * !n1->isBoundary;

	n2->cuy += mult2 * c23;
	n2->buy += mult2 * b23 * !n2->isBoundary;

	n4->auy += mult4 * c14 * !n4->isBoundary;
	n4->cuy += mult4 * b14;

	n3->auy += mult3 * c23 * !n3->isBoundary;
	n3->cuy += mult3 * b23;



	n1->cvy += (c14 + PRES_dtau / rho / dy) * mult1;
	n1->bvy += (b14 + PRES_dtau / rho / dy) * mult1 * !n1->isBoundary;
	n1->fvy += -mult1 * xp / rho * !n1->isBoundary;

	n2->cvy += (c23 + PRES_dtau / rho / dy) * mult2;
	n2->bvy += (b23 + PRES_dtau / rho / dy) * mult2 * !n2->isBoundary;
	n2->fvy += -mult2 * xp / rho * !n2->isBoundary;

	n4->avy += (c14 + PRES_dtau / rho / dy) * mult4 * !n4->isBoundary;
	n4->cvy += (b14 + PRES_dtau / rho / dy) * mult4;
	n4->fvy += mult4 * xp / rho * !n4->isBoundary;

	n3->avy += (c23 + PRES_dtau / rho / dy) * mult3 * !n3->isBoundary;
	n3->cvy += (b23 + PRES_dtau / rho / dy) * mult3;
	n3->fvy += mult3 * xp / rho * !n3->isBoundary;
}

void Cell::EvaluateXp(IterationStep itSt)
{
	if (itSt == xStep)
	{
		double hu14 = (n1->xu + n4->xu) / 2;
		double hu23 = (n3->xu + n2->xu) / 2;

		xp += -PRES_dtau * (hu23 - hu14) / dx * !isFreeOutlet;
	}
	else if (itSt == yStep)
	{
		double xv12 = (n1->xv + n2->xv) / 2;
		double xv34 = (n3->xv + n4->xv) / 2;

		xp += -PRES_dtau * (xv34 - xv12) / dy * !isFreeOutlet;
	}
}



void Cell::UpdatePressure(bool IsIterationsBreak, int n)
{
	if (IsIterationsBreak)
		p[n] = p_k;
	else
		p_k += xp;
}

void Cell::UpdateAverageVelocities()
{
	u12 = (n1->u_k + n2->u_k) / 2;
	u34 = (n3->u_k + n4->u_k) / 2;
	v14 = (n1->v_k + n4->v_k) / 2;
	v23 = (n2->v_k + n3->v_k) / 2;

	u14 = (n1->u_k + n4->u_k) / 2;
	u23 = (n2->u_k + n3->u_k) / 2;
	v12 = (n1->v_k + n2->v_k) / 2;
	v34 = (n3->v_k + n4->v_k) / 2;
}

double Cell::rot()
{
	double dvdx = (v23 - v14) / dx;
	double dudy = (u34 - u12) / dy;
	return dvdx - dudy;
}



