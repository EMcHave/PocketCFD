#include "pch.h"
#include "NavierStokes.h"
#include "NavierStokes.g.cpp"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <ppl.h>
#include <pplawait.h>

#define ALLNODES(m) for(int m = 0; m < NX * NY; m++)
#define ALLCELLS(m) for(int m = 0; m < (NX - 1) * (NY - 1); m++)

using namespace winrt::Windows::Foundation;

namespace winrt::CFD::implementation
{
	NavierStokes::NavierStokes()
		:NX(201), NY(101), l(1), h(0.5), dt(1.0/60), t(300.0/60), rho(1000), nu(pow(10, -6)), PRES_dtau(100)
	{
		this->dx = l / (NX - 1);
		this->dy = h / (NY - 1);
		this->dtau = dt;

		

		inlet_speed = 2.0;

		function<double(double, double)> BC1{
			[&](double y, double t) {return inlet_speed; }
		};
		function<double(double, double)> BC2{
		[&](double y, double t) {return inlet_speed; }
		};

		leftBC = BC1;
		rightBC = BC2;


	}

	NavierStokes::~NavierStokes()
	{
		if (cells != nullptr && nodes != nullptr)
			CleanSolution();
	}


	IAsyncActionWithProgress<Collections::IVector<double>> NavierStokes::Solve()
	{
		auto progress{ co_await winrt::get_progress_token() };
		winrt::apartment_context ui_thread;

		Collections::IVector<double> progressInfo(winrt::single_threaded_vector<double>());

		co_await concurrency::create_task([&] {SetInitialConditions(); });

		int it_counter = 0;

		double eps = pow(10, -2);		
		for (int n = 1; n < Nt(); n++)
		{

			SetBoundaryConditions(n);

			while (true)
			{
				co_await winrt::resume_background();

				ExplicitStep(n);

				MAX = MaxXi();
				bool IstimeToBreak = MAX < eps;



				if (IstimeToBreak)
				{
					ALLNODES(m) { Nod(m)->UpdateVelocities(IstimeToBreak, n); }
					ALLCELLS(m) { Cel(m)->UpdatePressure(IstimeToBreak); }
					break;
				}

				concurrency::parallel_for(0, NY - 1, 2, [&](int j)
					{
						for (int i = 0; i < NX - 1; i++)
						Cel(i, j)->XStepForCell(dt, dtau, rho, nu);
					});

				concurrency::parallel_for(1, NY - 1, 2, [&](int j)
					{
						for (int i = 0; i < NX - 1; i++)
						Cel(i, j)->XStepForCell(dt, dtau, rho, nu);
					});
				XStep(n);



				concurrency::parallel_for(0, NY - 1, 2, [&](int j)
					{
						for (int i = 0; i < NX - 1; i++)
						Cel(i, j)->YStepForCell(dt, dtau, rho, nu);
					});
				concurrency::parallel_for(1, NY - 1, 2, [&](int j)
					{
						for (int i = 0; i < NX - 1; i++)
						Cel(i, j)->YStepForCell(dt, dtau, rho, nu);
					});
				YStep();

				ALLCELLS(m) { Cel(m)->UpdatePressure(IstimeToBreak); }
				ALLNODES(m) { Nod(m)->UpdateVelocities(IstimeToBreak, n); }

				if (it_counter % 2 == 0)
				{
					co_await ui_thread;
					progressInfo.Append(n * 1.0 / (Nt() - 1));
					progressInfo.Append(MAX);
					progress(progressInfo);
					progressInfo.Clear();
				}
				it_counter++;
			}

			co_await ui_thread;
			progressInfo.Append(n * 1.0 / (Nt() - 1));
			progressInfo.Append(MAX);
			progress(progressInfo);
			progressInfo.Clear();
			
		}
	}

	void NavierStokes::TimeStep(int n, double eps)
	{
		
	}

	void NavierStokes::ExplicitStep(int n)
	{
		using namespace concurrency;


		ALLNODES(m)
		{
			Nod(m)->ResetNode(n);
		};

		parallel_for(0, NY - 1, 2, [&](int j)
			{
				for (int i = 0; i < NX - 1; i++)
				Cel(i, j)->EvaluateIntegralsForCell(dt, dtau, rho, nu, n);
			});

		parallel_for(1, NY - 1, 2, [&](int j)
			{
				for (int i = 0; i < NX - 1; i++)
				Cel(i, j)->EvaluateIntegralsForCell(dt, dtau, rho, nu, n);
			});

		ALLNODES(m)
		{
			Nod(m)->EvaluateXiuv(dtau, dt, n);
		};

	}


	void NavierStokes::XStep(int)
	{
		concurrency::parallel_for(1, NY - 1, [&](int j)
		{
			vector<double> au(NX), cu(NX), bu(NX), fu(NX);
			vector<double> av(NX), cv(NX), bv(NX), fv(NX);
			vector<double> hus, hvs;

			for (int i = 0; i < NX; i++)
			{
				au[i] = Nod(i, j)->aux;
				cu[i] = -Nod(i, j)->cux;
				bu[i] = Nod(i, j)->bux;
				fu[i] = -Nod(i, j)->fux;
				av[i] = Nod(i, j)->avx;
				cv[i] = -Nod(i, j)->cvx;
				bv[i] = Nod(i, j)->bvx;
				fv[i] = -Nod(i, j)->fvx;
			}
			hus = ThomasAlg(au, cu, bu, fu);
			hvs = ThomasAlg(av, cv, bv, fv);
			for (int i = 0; i < NX; i++)
			{
				Nod(i, j)->xu = hus[i];
				Nod(i, j)->xv = hvs[i];
			}
			for (int i = 0; i < NX; i++)
			{
				Nod(i, j)->fuy += hus[i] / dtau;
				Nod(i, j)->fvy += hvs[i] / dtau;
			}
		});
		ALLCELLS(m) { Cel(m)->EvaluateXp(IterationStep::xStep); }
	}



	void NavierStokes::YStep()
	{
		concurrency::parallel_for(1, NX - 1, [&](int i)
		{
			vector<double> au(NY), cu(NY), bu(NY), fu(NY);
			vector<double> av(NY), cv(NY), bv(NY), fv(NY);
			vector<double> xus, xvs;

			for (int j = NY - 1; j >= 0; j--)
			{
				au[j] = Nod(i, j)->auy;
				cu[j] = -Nod(i, j)->cuy;
				bu[j] = Nod(i, j)->buy;
				fu[j] = -Nod(i, j)->fuy;
				av[j] = Nod(i, j)->avy;
				cv[j] = -Nod(i, j)->cvy;
				bv[j] = Nod(i, j)->bvy;
				fv[j] = -Nod(i, j)->fvy;
			}
			xus = ThomasAlg(au, cu, bu, fu);
			xvs = ThomasAlg(av, cv, bv, fv);
			for (int j = NY - 1; j >= 0; j--)
			{
				Nod(i, j)->xu = xus[j];
				Nod(i, j)->xv = xvs[j];
			}
		});
		ALLCELLS(m) { Cel(m)->EvaluateXp(IterationStep::xStep); }
	}

	double NavierStokes::MaxXi()
	{
		double maxXi = 0;
		//#pragma omp parallel
		{
			//#pragma omp parallel for reduction(max:maxXi)
			ALLNODES(m)
			{
				maxXi = max(abs(Nod(m)->xu) / dtau, maxXi);
				maxXi = max(abs(Nod(m)->xv) / dtau, maxXi);
			}

			//#pragma omp parallel for reduction(max:maxXi)
			ALLCELLS(m) { maxXi = max(abs(Cel(m)->Divergence()), maxXi); }
		}

		return maxXi;
	}

	void NavierStokes::SetBoundaryConditions(int curMoment)
	{
		for (int j = 1; j < NY - 1; j++)
		{
			Nod(0, j)->isBoundary = true;
			Nod(NX - 1, j)->isBoundary = true;
		}
		for (int i = 0; i < NX; i++)
		{
			Nod(i, 0)->isBoundary = true;
			Nod(i, NY - 1)->isBoundary = true;
		}
		for (int j = 1; j < NY - 1; j++)
		{
			Nod(0, j)->u[curMoment] = leftBC(j * dy, curMoment * dt);
			Nod(NX - 1, j)->u[curMoment] = rightBC(j * dy, curMoment * dt);
			Nod(0, j)->u_k = Nod(0, j)->u[curMoment];
			Nod(NX - 1, j)->u_k = Nod(NX - 1, j)->u[curMoment];
		}

		for (int j = 45; j < 55; j++)
		{
			//Nod(10, j)->isBoundary = true;
			//Nod(11, j)->isBoundary = true;
			//Nod(12, j)->isBoundary = true;
			//Nod(8, j)->isBoundary = true;
			//Nod(9, j)->isBoundary = true;
		}
	}

	void NavierStokes::SetInitialConditions()
	{
		int Ncells = (NX - 1) * (NY - 1);
		int Nnodes = NX * NY;
		
		if (cells != nullptr && nodes != nullptr)
			CleanSolution();

		cells = new Cell * [Ncells];
		nodes = new Node * [Nnodes];

		int m = 0;
		for (int j = NY - 1; j >= 0; j--)
			for (int i = 0; i < NX; i++)
				*(nodes + m++) = new Node(i * dx, j * dy, m, this->Nt(), dt, dtau);
		double DX = dx;

		m = 0;
		for (int j = 0; j < NY - 1; j++)
			for (int i = 0; i < NX - 1; i++)
				*(cells + m++) = new Cell(
					Nod(i, j), Nod(i + 1, j),
					Nod(i + 1, j + 1), Nod(i, j + 1),
					m, this->Nt(), dx, dy, PRES_dtau);
		SetBoundaryConditions(0);
	}

	void NavierStokes::CleanSolution()
	{
		for (int i = 0; i < NX * NY; i++)
			delete nodes[i];
		for (int i = 0; i < (NX - 1) * (NY - 1); i++)
			delete cells[i];
		delete[] nodes;
		delete[] cells;
	}

	double NavierStokes::MaxVelocity(int n)
	{
		double max = 0;

		for (int m = 0; m < NX * NY; m++)
			if (nodes[m]->AbsVelocity(n) > max)
				max = nodes[m]->AbsVelocity(n);
		return max;
	}

	Node* const NavierStokes::Nod(int i, int j)
	{
		//assert(i < NX, j < NY);
		return nodes[j * NX + i];
	}

	Cell* const NavierStokes::Cel(int i, int j)
	{
		assert(i < NX - 1, j < NY - 1);
		return cells[j * (NX - 1) + i];
	}

	Node* const NavierStokes::Nod(int i)
	{
		assert(i < NX* NY);
		return nodes[i];
	}

	Cell* const NavierStokes::Cel(int i)
	{
		assert(i < (NX - 1)* (NY - 1));
		return cells[i];
	}

	event_token NavierStokes::PropertyChanged(Windows::UI::Xaml::Data::PropertyChangedEventHandler const& handler)
	{
		return m_propertyChanged.add(handler);
	}

	void NavierStokes::PropertyChanged(winrt::event_token const& token) noexcept
	{
		m_propertyChanged.remove(token);
	}
	
	event_token NavierStokes::IterationCompleted(Windows::Foundation::EventHandler<CFD::NavierStokesEventArgs> const& handler)
	{
		return m_iterationCompleted.add(handler);
	}

	void NavierStokes::IterationCompleted(winrt::event_token const& token) noexcept
	{
		m_iterationCompleted.remove(token);
	}
	
	double NavierStokes::Rho() { return rho; }

	void NavierStokes::Rho(double value)
	{
		rho = value;
	}

	double NavierStokes::Nu()
	{
		return nu;
	}

	void NavierStokes::Nu(double value)
	{
		nu = value;
		m_propertyChanged(*this, Windows::UI::Xaml::Data::PropertyChangedEventArgs{ L"Re" });
	}

	double NavierStokes::L()
	{
		return l;
	}

	void NavierStokes::L(double value)
	{
		l = value;
		dx = l / (NX - 1);
		m_propertyChanged(*this, Windows::UI::Xaml::Data::PropertyChangedEventArgs{ L"Re" });
	}

	double NavierStokes::H()
	{
		return h;
	}

	void NavierStokes::H(double value)
	{
		h = value;
		dy = h / (NY - 1);
	}

	int NavierStokes::Nx()
	{
		return NX;
	}

	void NavierStokes::Nx(int value)
	{
		NX = value;
		dx = l / (NX - 1);
	}

	int NavierStokes::Ny()
	{
		return NY;
	}

	void NavierStokes::Ny(int value)
	{
		NY = value;
		dy = h / (NY - 1);
	}

	double NavierStokes::T()
	{
		return t;
	}

	void NavierStokes::T(double  value)
	{
		t = value;
	}

	int NavierStokes::Nt()
	{
		return t / dt;
	}

	double NavierStokes::Pdtau()
	{
		return PRES_dtau;
	}

	void NavierStokes::Pdtau(double value)
	{
		PRES_dtau = value;
	}

	vector<double> NavierStokes::ThomasAlg(vector<double>& a, vector<double>& c,
		vector<double>& b, vector<double>& f)
	{
		int size = a.size();
		double* delta = new double[size];
		double* lambda = new double[size];
		vector<double> solution(size);


		delta[0] = -a[0] / c[0];
		lambda[0] = f[0] / c[0];
		double denom;
		int N1 = a.size() - 1;
		for (int i = 1; i < N1; i++)
		{
			denom = c[i] + a[i] * delta[i - 1];
			delta[i] = -b[i] / denom;
			lambda[i] = (f[i] - a[i] * lambda[i - 1]) / denom;
		}

		solution[N1] = (f[N1] - a[N1] * lambda[N1 - 1]) /
			(c[N1] + a[N1] * delta[N1 - 1]);

		for (int i = N1 - 1; i >= 0; i--)
			solution[i] = delta[i] * solution[i + 1] + lambda[i];

		delete[] delta;
		delete[] lambda;

		return solution;
	}


	void NavierStokes::XStep(vector<double>& au, vector<double>& cu,
		vector<double>& bu, vector<double>& fu,
		vector<double>& av, vector<double>& cv,
		vector<double>& bv, vector<double>& fv,
		vector<double>& hus, vector<double>& hvs, int n)
	{
		for (int j = 1; j < NY - 1; j++)
		{
			for (int i = 0; i < NX; i++)
			{
				au[i] = Nod(i, j)->aux;
				cu[i] = -Nod(i, j)->cux;
				bu[i] = Nod(i, j)->bux;
				fu[i] = -Nod(i, j)->fux;
				av[i] = Nod(i, j)->avx;
				cv[i] = -Nod(i, j)->cvx;
				bv[i] = Nod(i, j)->bvx;
				fv[i] = -Nod(i, j)->fvx;
			}
			hus = ThomasAlg(au, cu, bu, fu);
			hvs = ThomasAlg(av, cv, bv, fv);
			for (int i = 0; i < NX; i++)
			{
				Nod(i, j)->xu = hus[i];
				Nod(i, j)->xv = hvs[i];
			}
			for (int i = 0; i < NX; i++)
			{
				Nod(i, j)->fuy += hus[i] / dtau;
				Nod(i, j)->fvy += hvs[i] / dtau;
			}
		}
		ALLCELLS(m) { Cel(m)->EvaluateXp(IterationStep::xStep); }
	}

	void NavierStokes::YStep(vector<double>& au, vector<double>& cu,
		vector<double>& bu, vector<double>& fu,
		vector<double>& av, vector<double>& cv,
		vector<double>& bv, vector<double>& fv,
		vector<double>& xus, vector<double>& xvs)
	{
		for (int i = 1; i < NX - 1; i++)
		{
			for (int j = NY - 1; j >= 0; j--)
			{
				au[j] = Nod(i, j)->auy;
				cu[j] = -Nod(i, j)->cuy;
				bu[j] = Nod(i, j)->buy;
				fu[j] = -Nod(i, j)->fuy;
				av[j] = Nod(i, j)->avy;
				cv[j] = -Nod(i, j)->cvy;
				bv[j] = Nod(i, j)->bvy;
				fv[j] = -Nod(i, j)->fvy;
			}
			xus = ThomasAlg(au, cu, bu, fu);
			xvs = ThomasAlg(av, cv, bv, fv);
			for (int j = NY - 1; j >= 0; j--)
			{
				Nod(i, j)->xu = xus[j];
				Nod(i, j)->xv = xvs[j];
			}
		}
		ALLCELLS(m) { Cel(m)->EvaluateXp(IterationStep::yStep); }
	}
}
