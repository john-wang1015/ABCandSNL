#pragma once

#include "Params.hpp"
#include "Cell.hpp"
#include "delaunator.hpp"
#include <vector>
#include <assert.h>
#include <chrono>

using namespace::std;

double totalTime = 0;

class Pancreas
{
public:
	Pancreas(vector<Cell*> &input_cells, Params* parameters)
	{
		this->parameters = parameters;

		// Clone the input cells in case those cells are used in more than one particle
		map<Cell*, Cell*> map;

		for (Cell* cell : input_cells)
		{
			Cell* clone = new Cell(cell);
			cells.push_back(clone);
			map.insert(std::map<Cell*, Cell*>::value_type(cell, clone));
		}

		for (Cell* cell : cells)
			cell->updateSibling(map);

		src = dst = NULL;
	}

	~Pancreas()
	{
		//delete parameters; // is this right???

		for (Cell* cell : cells)
			delete cell;
	}

	Pancreas* CreateNewParticle(Params* parameters)
	{
		return new Pancreas(cells, parameters);
	}

private:
	vector<Cell*> cells;
	Params* parameters;
	vector<Cell*> new_cells; // accumulated during the latest iteration
	vector<Cell*> boundaryCells;
	Cell* src, *dst; // the two cells that are furthest apart - using for computing dimensions of tumour

	void DetermineNeighbours()
	{
		double* coords = new double[cells.size() * 2];
		int i = 0;
		for (Cell* cell : cells)
		{
			cell->clearNeighbours();
			coords[i++] = cell->currentState.X;
			coords[i++] = cell->currentState.Y;
		}

		Delaunator TRI(coords, (int)cells.size());

		delete[] coords;

        #define nextHalfedge(e) ((e % 3 == 2) ? e - 2 : e + 1)

		for (size_t e = 0; e < TRI.trianglesLen; e++)
			if (TRI.halfedges[e] == -1 || e > TRI.halfedges[e])
			{
				Cell* P = cells[TRI.triangles[e]];
				Cell* Q = cells[TRI.triangles[nextHalfedge(e)]];
				P->Neighbours.push_back(Q);
				Q->Neighbours.push_back(P);
			}
	}

	bool HealthyCellsBeyondRadius(double radius)
	{
		for (Cell* cell : cells)
			if (cell->currentState.type == CellType::Healthy && cell->DistanceSquaredFromCentre() >= Sqr(radius))
				return true;
		return false;
	}

	void AddNewCell(Cell* new_cell)
	{
		if (new_cell != NULL)
			new_cells.push_back(new_cell);
	}

	void AddMoreTissue(double moving_rim, double max_tumour_radius)
	{
		const double Y_spacing = 3;
		double X_spacing = Y_spacing * sqrt(3) / 2;

		double new_radius = max_tumour_radius + moving_rim + 20;
		int i = 0;
		for (double X = -new_radius; X <= new_radius; X += X_spacing)
		{
			// should offset odd columns by 1 or Y_spacing/2 ???
			for (double Y = -new_radius + (i % 2) * Y_spacing / 2; Y <= new_radius; Y += Y_spacing)
			{
				double distance_squared = DistanceSquared(X, Y, 0, 0);
				if (distance_squared > Sqr(max_tumour_radius) && distance_squared <= Sqr(max_tumour_radius + moving_rim + 10)) // why +10 vs +20 above?
					AddNewCell(new Cell(X, Y, Params::s, NULL, CellType::Healthy, parameters->gage));
			}
			i++;
		}
	}
	void MoreTissueAddedIfNecessary()
	{
		const int moving_rim = 10;
		double tumour_radius = TumourRadius();
		if (!HealthyCellsBeyondRadius(tumour_radius + moving_rim))
			AddMoreTissue(moving_rim, tumour_radius);
	}

	void DetermineBoundaryCells()
	{
		if (cells[0]->Neighbours.size() == 0)
			DetermineNeighbours();

		boundaryCells.clear();
		for (Cell* cell : cells)
			if (cell->currentState.type != CellType::Healthy)
			{
				for (Cell* neighbour: cell->Neighbours)
					if (neighbour->currentState.type == CellType::Healthy)
					{
						boundaryCells.push_back(cell);
						break;
					}
			}
	}

	// Still used for determining when to add more cells
	double TumourRadius()
	{
		double maxDistance = 0;
		for (Cell* cell : cells)
			if (cell->currentState.type != CellType::Healthy)
			{
				double distance = cell->DistanceSquaredFromCentre();
				if (distance > maxDistance)
					maxDistance = distance;
			}
		return sqrt(maxDistance);
	}

	double DistanceToLine(Cell* cell)
	{
		return (dst->currentState.X - src->currentState.X)  * (src->currentState.Y - cell->currentState.Y) - 
			   (src->currentState.X - cell->currentState.X) * (dst->currentState.Y - src->currentState.Y);
	}

	void SimulateOneHour()
	{
		DetermineBoundaryCells();

		for (Cell* cell : cells)
		{
			cell->Renew(); // set new state to current state

			if (cell->OnBoundary())
				cell->PossiblyPSCInfectNeighbour(parameters);

			if (cell->currentState.type == CellType::Healthy)
				cell->Move();
			else
			{
				if (cell->currentState.age < parameters->gage)
					cell->LengthenSpring(parameters);

				// cancer cells either proliferate or move ...
				Cell* newCell = cell->PossiblyPoliferate(boundaryCells, parameters);
				if (newCell != NULL)
					AddNewCell(newCell);
				else
					cell->Move();
			}
		}

		for (Cell* cell : cells)
			cell->UpdateState(); // set current state to new state

		MoreTissueAddedIfNecessary();

		// Only add new cells after all processing for this hour is complete.
		cells.insert(cells.end(), new_cells.begin(), new_cells.end());
		new_cells.clear();

		DetermineNeighbours();
	}

public:
	double TumourVolume()
	{
		if (boundaryCells.size() == 0)
			DetermineBoundaryCells();

		if (boundaryCells.size() <= 1)
			return 0;

		double longest = -1;

		// Find the two cells that are furthest apart from each other
		for (Cell* cell1 : boundaryCells)
		{
			for (Cell* cell2 : boundaryCells)
			{
				double distSquared = cell1->DistanceSquaredTo(cell2);
				if (distSquared > longest)
				{
					longest = distSquared;
					src = cell1;
					dst = cell2;
				}
			}
		}

		double length = sqrt(longest);

		// Find width from cells furthest from the centre line
		double minDist = 0, maxDist = 0;
		for (Cell* cell : boundaryCells)
		{
			double distance = DistanceToLine(cell);
			if (distance < minDist) minDist = distance;
			if (distance > maxDist) maxDist = distance;
		}

		double scale = sqrt(Sqr(dst->currentState.X - src->currentState.X) + Sqr(dst->currentState.Y - src->currentState.Y));
		minDist /= scale;
		maxDist /= scale;

		double width = maxDist - minDist;

		assert(width <= length);

		width *= 0.1728;
		length *= 0.1728;

		double volume = width * width * length / 2;
		//assert(!isnan(volume));
		return volume;
	}

	void CreateInitialTumour()
	{
		AddMoreTissue(10, 0);
		// append the new cells to the current cells ...
		cells.insert(cells.end(), new_cells.begin(), new_cells.end());
		new_cells.clear();

		Cell* closestDistanceToCentre = NULL;
		double closestDistance = MAX_DBL;
		for (Cell* cell : cells)
		{
			double distance = cell->DistanceSquaredFromCentre();
			if (distance < closestDistance)
			{
				closestDistanceToCentre = cell;
				closestDistance = distance;
			}
		}
		closestDistanceToCentre->Infect();
	}

	double SimulateOneDay(int day)
	{
		DetermineNeighbours();

		for (int hour = 1; hour <= Params::tinterval; hour++)
		{
			//parameters->created = parameters->opportunities = 0;
			SimulateOneHour();
			//printf("day %d, hour %d, %d of %d created\n", day, hour, parameters->created, parameters->opportunities);
		}
		return TumourVolume();
	}
	void UpdateParameters(Params* parameters)
    {
        this->parameters = parameters;
    }
};
