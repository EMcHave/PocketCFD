#pragma once
#include <functional>
#include "NavierStokes.g.h"
#include "Cell.h"
#include "NavierStokesEventArgs.g.h"

using namespace std;
using namespace winrt::Windows::Foundation;
using namespace winrt::Windows::Foundation::Collections;

namespace winrt::CFD::implementation
{
	enum class TypeOfBC
	{
		Wall,
		Speed,
		Pressure
	};
	enum class Boundary
	{
		Left,
		Right, 
		Top, 
		Bottom
	};
	struct BoundaryCondition
	{
		TypeOfBC Type;
		double Vx;
		double Vy;
		double P;
	};
	struct NavierStokes : NavierStokesT<NavierStokes>
	{

		//Cell** cells;
		//Node** nodes;

		vector<Cell*> cells;
		vector<Node*> nodes;

		NavierStokes();
		~NavierStokes();


		IAsyncActionWithProgress<IVector<double>> Solve(map<Boundary, BoundaryCondition>, vector<Point>);
		void SetBoundaryConditions(int);
		void SetInitialConditions(vector<Point>);
		void CleanSolution();
		double MaxVelocity(int n);
		double MaxPressure(int n);
		double MinPressure(int n);
;		Node* const Nod(int i, int j);
		Cell* const Cel(int i, int j);
		Node* const Nod(int i);
		Cell* const Cel(int i);

		bool Solved() { return solved; }
		void Solved(bool s) { solved = s; }

		event_token PropertyChanged(Windows::UI::Xaml::Data::PropertyChangedEventHandler const& handler);
		void PropertyChanged(winrt::event_token const& token) noexcept;
		event_token IterationCompleted(EventHandler<CFD::NavierStokesEventArgs> const& handler);
		void IterationCompleted(winrt::event_token const& token) noexcept;

		double Rho() { return rho; }
		void Rho(double value) { rho = value; }

		double Nu() { return nu; }
		void Nu(double);

		double L() { return l; }
		void L(double);

		double H() { return h; }
		void H(double);

		int Nx() { return NX; }
		void Nx(int);

		int Ny() { return NY; }
		void Ny(int);

		double Dx() { return L() / (Nx() - 1); }
		double Dy() { return H() / (Ny() - 1); }

		double T() { return t; }
		void T(double value) { t = value; }

		int Nt() { return t / dt; }

		double Dt() { return dt; }
		void Dt(double value) { dt = value; }

		double Pdtau() { return PRES_dtau; }
		void Pdtau(double value) { PRES_dtau = value; }

		double Dtau() { return dtau; }
		void Dtau(double value) { dtau = value; }

		double Re() { return inlet_speed * L() / Nu(); }

		double Eps() { return eps; }
		void Eps(double value) { eps = value; }

		void ThomasAlg(vector<double>& a, vector<double>& с,
			vector<double>& b, vector<double>& f, vector<double>& solution);
	private:
		double rho;
		double nu;

		double l;
		double h;
		double dx;
		double dy;
		int NX;
		int NY;

		double t;
		int NT;
		double dt;
		double dtau;
		double PRES_dtau;
		double eps;
		double inlet_speed;

		double MAX;
		double maxXiU, maxXiV, maxXiP;

		bool solved;

		map<Boundary, BoundaryCondition> BCs;
		
		event<Windows::UI::Xaml::Data::PropertyChangedEventHandler> m_propertyChanged;
		event<EventHandler<CFD::NavierStokesEventArgs>> m_iterationCompleted;

		//function<double(double, double)> leftBC;
		//function<double(double, double)> rightBC;

		void ExplicitStep(int n);
		void XStep(int);
		void YStep();

		double MaxXi();
		double MaxXiU();
		double MaxXiV();
		double MaxXiP();
	};


	struct NavierStokesEventArgs : NavierStokesEventArgsT<NavierStokesEventArgs>
	{
		NavierStokesEventArgs() = default;
		NavierStokesEventArgs(double max_Xi):m_maximumXi(max_Xi)
		{}
		double MaximumXi() { return m_maximumXi; }
	private:
		double m_maximumXi{ 0.0 };
	};

}
namespace winrt::CFD::factory_implementation
{
	struct NavierStokes : NavierStokesT<NavierStokes, implementation::NavierStokes>
	{
	};
}
