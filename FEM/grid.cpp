#include "grid.h"
#include <vector>
#include <fstream>
Grid::Grid()
{
	std::ifstream file("grid.txt", std::ios::in);
	if (!file.good())
	{
		throw std::runtime_error("Error with opening file. Check path or filename.");
	}
	file >> H;
	file >> L;
	file >> nH;
	file >> nL;
	file >> t0;
	nN = nH*nL;
	nE = (nH - 1)*(nL - 1);
	nodes.reserve(nN);
	elements.reserve(nE);
}

Grid::~Grid()
{
}

void Grid::setGrid()
{
	float x = 0.0, y = 0.0, dx=L/nL, dy=H/nH;
	bool L=true;
	for (size_t i = 1; i <= nN; i++)
	{
		Node mynode(i, x, y, t0);
		nodes.push_back(mynode);
		y += dy;
		if (y == H)
		{
			y = 0.0;
			x += dx;
		}
		if (i <= nN - nH)
		{
			if (i%nH != 0)
			{
				std::array<unsigned int, 4> IDs;
				IDs[0] = i;
				IDs[1] = i + 1;
				IDs[2] = i + nH;
				IDs[3] = i + nH + 1;
				Element myNewElement(IDs);
				elements.push_back(myNewElement);
			}
		}
	}
}
void Grid::printGrid()
{
	std::cout << "Nodes:" << "\n";
	for (auto&x:nodes)
	{
		std::cout << x << "\n";
	}
	std::cout << "Elements:" << "\n";
	for (auto&x:elements)
	{
		for (auto&y:x.getArray())
		{
			std::cout << y << " ";
		}
		std::cout << "\n";
	}
}
void Grid::CalculateH()
{
	std::array<double, 4> cordsx;
	std::array<double, 4> cordsy;
	for (size_t i = 0; i < 4; ++i) // we need to get cords (x,y) for our element
	{
		unsigned int tmp = elements[0].getID(i)-1; // -1 becaouse ID-s are 1-4 no 0-3 :(
		cordsx[i] = nodes[tmp].getX();
		cordsy[i] = nodes[tmp].getY();

	}
	this->elements[0].calculH(cordsx,cordsy);
}