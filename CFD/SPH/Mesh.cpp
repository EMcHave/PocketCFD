#include "Mesh.h"

Ice::Mesh::Mesh(float2 size, float radius) : m_size(size), m_cellSize(radius)
{
	float dR = radius;
	NX = (int)size.x / dR;
	NY = (int)size.y / dR;

	float2 sdvig(size.x * 0.5, size.y * 0.5);

	for (int j = 0; j < NY; j++)
		for (int i = 0; i < NX; i++)
				cells.push_back(std::make_shared<Ice::Cell>(i, j, dR, sdvig));
}

void Ice::Mesh::CleanCells()
{
	for (auto cell : cells)
		cell->Clear();
}

std::shared_ptr<Ice::Cell> Ice::Mesh::CellByPosition(float2 pos)
{
	unsigned int i = (pos.x + 0.5 * m_size.x) / m_cellSize;
	unsigned int j = (pos.y + 0.5 * m_size.y) / m_cellSize;

	if (i < NX && j < NY)
		return Cel(i, j);
	else
		return std::shared_ptr<Ice::Cell>(nullptr);
}

std::shared_ptr<Ice::Cell> Ice::Mesh::Cel(int i, int j, int k)
{
	return cells.at(NX * NZ * j + NZ * i + k);
}

Ice::Cell::Cell(int i, int j, int k, float r, float2 sdvig)
{
	isEmpty = true;
}

void Ice::Cell::AddParticle(SParticle p)
{
	particles.push_back(p);
	if(isEmpty)
		isEmpty = false;
}

void Ice::Cell::Clear()
{
	particles.clear();
	isEmpty = true;
}
