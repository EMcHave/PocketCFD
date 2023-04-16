#include "pch.h"
#include "NavierStokes.h"
#include "NavierStokes.g.cpp"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <ppl.h>
#include <pplawait.h>
#include <type_traits>

#define ALLNODES(m) for(int m = 0; m < NX * NY; m++)
#define ALLCELLS(m) for(int m = 0; m < (NX - 1) * (NY - 1); m++)

using namespace winrt::Windows::Foundation;

namespace winrt::CFD::implementation
{
	NavierStokes::NavierStokes()
		:NX(101), NY(51), l(1), h(0.5), dt(1.0/60), t(300.0/60), rho(1000), nu(pow(10, -6)), PRES_dtau(1000)
	{
		this->dx = l / (NX - 1);
		this->dy = h / (NY - 1);
		this->dtau = dt;
		this->eps = pow(10, -2);
		

		inlet_speed = 1.0;

		function<double(double, double)> BC1{
			[&](double y, double t) {return inlet_speed; }
		};
		function<double(double, double)> BC2{
		[&](double y, double t) {return inlet_speed; }
		};
	}

	NavierStokes::~NavierStokes()
	{
		if (!cells.empty() && !nodes.empty())
			CleanSolution();
	}


	IAsyncActionWithProgress<Collections::IVector<double>> NavierStokes::Solve(map<Boundary, BoundaryCondition> BCs, vector<Point> boundaries)
	{
		this->BCs = BCs;

		auto progress{ co_await winrt::get_progress_token() };
		winrt::apartment_context ui_thread;

		Collections::IVector<double> progressInfo(winrt::single_threaded_vector<double>());

		co_await concurrency::create_task([&] {SetInitialConditions(boundaries); });


		int it_counter = 0;

		for (int n = 1; n < Nt(); n++)
		{
			SetBoundaryConditions(n);

			while (true)
			{
				co_await winrt::resume_background();

				ExplicitStep(n);
				maxXiU = MaxXiU();
				maxXiV = MaxXiV();
				maxXiP = MaxXiP();
				MAX = max(maxXiU, max(maxXiV, maxXiP));
				bool IstimeToBreak = MAX < eps;



				if (IstimeToBreak)
				{
					ALLNODES(m) { Nod(m)->UpdateVelocities(IstimeToBreak, n); }
					ALLCELLS(m) { Cel(m)->UpdatePressure(IstimeToBreak, n); }
					break;
				}


				///////////////// X - STEP ///////////////////

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


				if (BCs[Boundary::Right].Type == TypeOfBC::Pressure)
					for (int j = 0; j < NY; j++)
					{
						Nod(NX - 1, j)->aux = Nod(NX - 1, j)->cux;
						Nod(NX - 1, j)->avx = Nod(NX - 1, j)->cvx;
					}
				if (BCs[Boundary::Left].Type == TypeOfBC::Pressure)
					for (int j = 0; j < NY; j++)
					{
						Nod(0, j)->bux = Nod(0, j)->cux;
						Nod(0, j)->bvx = Nod(0, j)->cvx;
					}
				XStep(n);


				if (BCs[Boundary::Top].Type == TypeOfBC::Pressure)
					for (int i = 0; i < NX; i++) {
						Nod(i, NY - 1)->xu = Nod(i, NY - 2)->xu;
						Nod(i, NY - 1)->xv = Nod(i, NY - 2)->xv;
					}

				if (BCs[Boundary::Bottom].Type == TypeOfBC::Pressure)
					for (int i = 0; i < NX; i++) {
						Nod(i, 0)->xu = Nod(i, 1)->xu;
						Nod(i, 0)->xv = Nod(i, 1)->xv;
					}







				///////////////// Y - STEP ///////////////////

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
				
				
				if (BCs[Boundary::Top].Type == TypeOfBC::Pressure)
					for (int i = 0; i < NX; i++)
					{
						Nod(i, NY - 1)->auy = Nod(i, NY - 1)->cuy;
						Nod(i, NY - 1)->avy = Nod(i, NY - 1)->cvy;
						
					}
				if (BCs[Boundary::Bottom].Type == TypeOfBC::Pressure)
					for (int i = 0; i < NX; i++)
					{
						Nod(i, 0)->buy = Nod(i, 0)->cuy;
						Nod(i, 0)->bvy = Nod(i, 0)->cvy;
					}

				YStep();

				
				if (BCs[Boundary::Right].Type == TypeOfBC::Pressure)
					for (int j = 0; j < NY; j++) {
						Nod(NX - 1, j)->xu = Nod(NX - 2, j)->xu;
						Nod(NX - 1, j)->xv = Nod(NX - 2, j)->xv;
					}

				if (BCs[Boundary::Left].Type == TypeOfBC::Pressure)
					for (int j = 0; j < NY; j++) {
						Nod(0, j)->xu = Nod(1, j)->xu;
						Nod(0, j)->xv = Nod(1, j)->xv;
					}

				ALLNODES(m) { Nod(m)->UpdateVelocities(IstimeToBreak, n); }
				ALLCELLS(m) { Cel(m)->UpdatePressure(IstimeToBreak, n); }

				if (it_counter % 20 == 0)
				{
					co_await ui_thread;
					progressInfo.Append(n * 1.0 / (Nt() - 1));
					progressInfo.Append(maxXiU);
					progressInfo.Append(maxXiV);
					progressInfo.Append(maxXiP);
					progressInfo.Append(MAX);
					progress(progressInfo);
					progressInfo.Clear();
				}
				it_counter++;
			}
			co_await ui_thread;
			progressInfo.Append(n * 1.0 / (Nt() - 1));
			progressInfo.Append(maxXiU);
			progressInfo.Append(maxXiV);
			progressInfo.Append(maxXiP);
			progressInfo.Append(MAX);
			progress(progressInfo);
			progressInfo.Clear();
			
		}
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
			vector<double> hus(NX), hvs(NX);

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
			ThomasAlg(au, cu, bu, fu, hus);
			ThomasAlg(av, cv, bv, fv, hvs);
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
			vector<double> xus(NY), xvs(NY);

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
			ThomasAlg(au, cu, bu, fu, xus);
			ThomasAlg(av, cv, bv, fv, xvs);
			for (int j = NY - 1; j >= 0; j--)
			{
				Nod(i, j)->xu = xus[j];
				Nod(i, j)->xv = xvs[j];
			}
		});
		ALLCELLS(m) { Cel(m)->EvaluateXp(IterationStep::yStep); }
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

	double NavierStokes::MaxXiU()
	{
		double maxXi = 0;
		ALLNODES(m) { maxXi = max(abs(Nod(m)->xu) / dtau, maxXi); }
		return maxXi;
	}

	double NavierStokes::MaxXiV()
	{
		double maxXi = 0;
		ALLNODES(m) { maxXi = max(abs(Nod(m)->xv) / dtau, maxXi); }
		return maxXi;
	}

	double NavierStokes::MaxXiP()
	{
		double maxXi = 0;
		ALLCELLS(m) { maxXi = max(abs(Cel(m)->Divergence()), maxXi); }
		return maxXi;
	}

	void NavierStokes::SetBoundaryConditions(int curMoment)
	{

		switch (BCs[Boundary::Left].Type)
		{
		case TypeOfBC::Speed:
		{
			for (int j = 1; j < NY - 1; j++)
			{
				Nod(0, j)->u[curMoment] = BCs[Boundary::Left].Vx;
				Nod(0, j)->v[curMoment] = BCs[Boundary::Left].Vy;
			}
			break;
		}
		case TypeOfBC::Pressure:
		{
			for (int j = 0; j < NY - 1; j++)
			{
				Cel(0, j)->p[curMoment] = BCs[Boundary::Left].P;
				Cel(0, j)->p_k = Cel(0, j)->p[curMoment];
			}
			break;
		}
		default:
			break;
		}

		switch (BCs[Boundary::Right].Type)
		{
		case TypeOfBC::Speed:
		{
			for (int j = 1; j < NY - 1; j++)
			{
				Nod(NX - 1, j)->u[curMoment] = BCs[Boundary::Right].Vx;
				Nod(NX - 1, j)->v[curMoment] = BCs[Boundary::Right].Vy;
			}
			break;
		}
		case TypeOfBC::Pressure:
		{
			for (int j = 0; j < NY - 1; j++)
				Cel(NX - 2, j)->p[curMoment] = BCs[Boundary::Right].P;
			break;
		}
		default:
			break;
		}

		switch (BCs[Boundary::Top].Type)
		{
		case TypeOfBC::Speed:
		{
			for (int i = 1; i < NX - 1; i++)
			{
				Nod(i, NY - 1)->u[curMoment] = BCs[Boundary::Top].Vx;
				Nod(i, NY - 1)->v[curMoment] = BCs[Boundary::Top].Vy;
			}
			break;
		}
		case TypeOfBC::Pressure:
		{
			for (int i = 0; i < NX - 1; i++)
				Cel(i, NY - 2)->p[curMoment] = BCs[Boundary::Top].P;
			break;
		}
		default:
			break;
		}

		switch (BCs[Boundary::Bottom].Type)
		{
		case TypeOfBC::Speed:
		{
			for (int i = 1; i < NX - 1; i++)
			{
				Nod(i, 0)->u[curMoment] = BCs[Boundary::Bottom].Vx;
				Nod(i, 0)->v[curMoment] = BCs[Boundary::Bottom].Vy;
			}
			break;
		}
		case TypeOfBC::Pressure:
		{
			for (int i = 0; i < NX - 1; i++)
				Cel(i, 0)->p[curMoment] = BCs[Boundary::Bottom].P;
			break;
		}
		default:
			break;
		}
	}

	void NavierStokes::SetInitialConditions(vector<Point> boundaries)
	{		
		if (!cells.empty() && !nodes.empty())
			CleanSolution();

		cells = vector<Cell*>((NX - 1) * (NY - 1), nullptr);
		nodes = vector<Node*>(NX * NY, nullptr);

		int m = 0;
		for (int j = NY - 1; j >= 0; j--)
			for (int i = 0; i < NX; i++)
				nodes[m++] = new Node(i * dx, j * dy, m, this->Nt(), dt, dtau);

		m = 0;
		for (int j = 0; j < NY - 1; j++)
			for (int i = 0; i < NX - 1; i++)
				cells[m++] = new Cell(
					Nod(i, j), Nod(i + 1, j),
					Nod(i + 1, j + 1), Nod(i, j + 1),
					m, this->Nt(), dx, dy, PRES_dtau);
		SetBoundaryConditions(0);

		/////////// Границы сверху и снизу ///////////
		for (int i = 0; i < NX; i++)
		{
			Nod(i, 0)->isBoundary = true;
			Nod(i, 0)->v_k = Nod(i, 0)->v[0];
			Nod(i, 0)->u_k = Nod(i, 0)->u[0];
			Nod(i, NY - 1)->isBoundary = true;
			Nod(i, NY - 1)->v_k = Nod(i, NY - 1)->v[0];
			Nod(i, NY - 1)->u_k = Nod(i, NY - 1)->u[0];
		}
		/////////// Стенки слева и справа ////////////
		for (int j = 0; j < NY; j++)
		{
			Nod(0, j)->isBoundary = true;
			Nod(0, j)->u_k = Nod(0, j)->u[0];
			Nod(0, j)->v_k = Nod(0, j)->v[0];
			Nod(NX - 1, j)->isBoundary = true;
			Nod(NX - 1, j)->u_k = Nod(NX - 1, j)->u[0];
			Nod(NX - 1, j)->v_k = Nod(NX - 1, j)->v[0];
		}

		/////////// Давление ///////////////
		
		if(BCs[Boundary::Right].Type == TypeOfBC::Pressure)
			for (int j = 0; j < NY - 1; j++)
			{
				Cel(NX - 2, j)->isFreeOutlet = true;
				Cel(NX - 2, j)->p_k = Cel(NX - 2, j)->p[0];
			}
		if (BCs[Boundary::Left].Type == TypeOfBC::Pressure)
			for (int j = 0; j < NY - 1; j++)
			{
				Cel(0, j)->isFreeOutlet = true;
				Cel(0, j)->p_k = Cel(0, j)->p[0];
			}
		if (BCs[Boundary::Top].Type == TypeOfBC::Pressure)
			for (int i = 0; i < NX - 1; i++)
			{
				Cel(i, NY - 2)->isFreeOutlet = true;
				Cel(i, NY - 2)->p_k = Cel(i, NY - 2)->p[0];
			}
		if (BCs[Boundary::Bottom].Type == TypeOfBC::Pressure)
			for (int i = 0; i < NX - 1; i++)
			{
				Cel(i, 0)->isFreeOutlet = true;
				Cel(i, 0)->p_k = Cel(i, 0)->p[0];
			}


		for (auto node : boundaries)
			Nod((int)node.X, (int)node.Y)->isBoundary = true;
	}

	void NavierStokes::CleanSolution()
	{
		/*
		for (int i = 0; i < NX * NY; i++)
			if(dynamic_cast<Node*>(nodes[i]))
				delete nodes[i];
		for (int i = 0; i < (NX - 1) * (NY - 1); i++)
			if(cells[i] != NULL)
				delete cells[i];
		delete[] nodes;
		delete[] cells;

		nodes = nullptr;
		cells = nullptr;
		*/

		for (int i = 0; i < nodes.size(); i++)
			delete nodes[i];
		for (int i = 0; i < cells.size(); i++)
			delete cells[i];

		cells.clear();
		nodes.clear();
	}

	double NavierStokes::MaxVelocity(int n)
	{
		double max = 0;

		for (int m = 0; m < NX * NY; m++)
			if (nodes[m]->AbsVelocity(n) > max)
				max = nodes[m]->AbsVelocity(n);
		return max;
	}

	double NavierStokes::MaxPressure(int n)
	{
		double max = 0;

		for (int m = 0; m < (NX - 1) * (NY - 1); m++)
			if (cells[m]->p[n] > max)
				max = cells[m]->p[n];
		return max;
	}

	double NavierStokes::MinPressure(int n)
	{
		double min = cells[0]->p[n];

		for (int m = 0; m < (NX - 1) * (NY - 1); m++)
			if (cells[m]->p[n] < min)
				min = cells[m]->p[n];
		return min;
	}

	Node* const NavierStokes::Nod(int i, int j)
	{
		assert(i < NX, j < NY);
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

	void NavierStokes::Nu(double value)
	{
		nu = value;
		m_propertyChanged(*this, Windows::UI::Xaml::Data::PropertyChangedEventArgs{ L"Re" });
	}

	void NavierStokes::L(double value)
	{
		l = value;
		dx = l / (NX - 1);
		m_propertyChanged(*this, Windows::UI::Xaml::Data::PropertyChangedEventArgs{ L"Re" });
	}

	void NavierStokes::H(double value)
	{
		h = value;
		dy = h / (NY - 1);
	}

	void NavierStokes::Nx(int value)
	{
		NX = value;
		dx = l / (NX - 1);
	}

	void NavierStokes::Ny(int value)
	{
		NY = value;
		dy = h / (NY - 1);
	}

	void NavierStokes::ThomasAlg(vector<double>& a, vector<double>& c,
		vector<double>& b, vector<double>& f, vector<double>& solution)
	{
		int size = a.size();
		double* delta = new double[size];
		double* lambda = new double[size];
		


		delta[0] = -b[0] / c[0];
		lambda[0] = f[0] / c[0];
		double denom;
		
		int N = a.size();
		int N1 = N - 1;
		for (int i = 1; i < N1; i++)
		{
			denom = c[i] + a[i] * delta[i - 1];
			delta[i] = -b[i] / denom;
			lambda[i] = (f[i] - a[i] * lambda[i - 1]) / denom;
		}

		solution[N1] = (f[N1] - a[N1] * lambda[N1 - 1]) /
			(c[N1] + a[N1] * delta[N1 - 1]);

		//solution[N1] = lambda[N1];

		for (int i = N1 - 1; i >= 0; i--)
			solution[i] = delta[i] * solution[i + 1] + lambda[i];
		

		delete[] delta;
		delete[] lambda;
	}
}
