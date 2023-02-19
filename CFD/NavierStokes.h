#pragma once
#include <functional>
#include "NavierStokes.g.h"
#include "Cell.h"
#include "NavierStokesEventArgs.g.h"

using namespace std;
using namespace winrt::Windows::Foundation;

namespace winrt::CFD::implementation
{
	struct NavierStokes : NavierStokesT<NavierStokes>
	{

		Cell** cells;
		Node** nodes;

		NavierStokes();
		~NavierStokes();


		IAsyncActionWithProgress<Collections::IVector<double>> Solve();
		void SetBoundaryConditions(int);
		void SetInitialConditions();
		void CleanSolution();
		double MaxVelocity(int n);
		Node* const Nod(int i, int j);
		Cell* const Cel(int i, int j);
		Node* const Nod(int i);
		Cell* const Cel(int i);

		bool Solved() { return solved; }
		void Solved(bool s) { solved = s; }

		event_token PropertyChanged(Windows::UI::Xaml::Data::PropertyChangedEventHandler const& handler);
		void PropertyChanged(winrt::event_token const& token) noexcept;
		event_token IterationCompleted(EventHandler<CFD::NavierStokesEventArgs> const& handler);
		void IterationCompleted(winrt::event_token const& token) noexcept;

		double Rho();
		void Rho(double);

		double Nu();
		void Nu(double);

		double L();
		void L(double);

		double H();
		void H(double);

		int Nx();
		void Nx(int);

		int Ny();
		void Ny(int);

		double Dx() { return L() / (Nx() - 1); }
		double Dy() { return H() / (Ny() - 1); }

		double T();
		void T(double);

		int Nt();

		double Dt() { return dt; }
		void Dt(double value) { dt = value; }

		double Pdtau();
		void Pdtau(double);

		double Re() { return inlet_speed * L() / Nu(); }

		vector<double> ThomasAlg(vector<double>& a, vector<double>& с,
			vector<double>& b, vector<double>& f);
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
		double inlet_speed;

		double MAX;

		bool solved;
		
		event<Windows::UI::Xaml::Data::PropertyChangedEventHandler> m_propertyChanged;
		event<EventHandler<CFD::NavierStokesEventArgs>> m_iterationCompleted;

		function<double(double, double)> leftBC;
		function<double(double, double)> rightBC;

		void TimeStep(int n, double eps);
		void ExplicitStep(int n);
		void XStep(vector<double>&, vector<double>&, vector<double>&, vector<double>&,
			vector<double>&, vector<double>&, vector<double>&, vector<double>&,
			vector<double>&, vector<double>&, int);
		void XStep(int);
		void YStep(vector<double>&, vector<double>&, vector<double>&, vector<double>&,
			vector<double>&, vector<double>&, vector<double>&, vector<double>&,
			vector<double>&, vector<double>&);
		void YStep();

		double MaxXi();
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
