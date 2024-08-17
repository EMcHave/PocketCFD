#pragma once
#include "SParticle.h"
#include <vector>
#include <memory>

namespace Ice
{
	struct Cell
	{
		Cell(int i, int j, float r, float2 sdvig);
		std::vector<SParticle>& Particles() { return particles; }
		void AddParticle(Particle p);
		void Clear();
		bool isEmpty;
	private:
		std::vector<SParticle> particles;
	};

	class Mesh
	{
		std::vector<std::shared_ptr<Ice::Cell>> cells;
		float2 m_size;

		float m_cellSize;
		int NX;
		int NY;
		int NZ;
	public:
		Mesh(float2 size, float radius);
		
		void CleanCells();
		float2 Size() { return float2(NX, NY); }
		std::vector<std::shared_ptr<Ice::Cell>>& Cells() { return cells; }
		
		std::shared_ptr<Ice::Cell> CellByPosition(float2 pos);
		std::shared_ptr<Ice::Cell> Cel(int i, int j, int k);
	};
}


