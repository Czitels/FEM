#include "grid.h"
#include "element.h"
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
	float x = 0.0, y = 0.0, dx=L/(nL-1), dy=H/(nH-1);
	for (size_t i = 1; i <= nN; i++)
	{
		Node mynode(i, x, y, t0); // -1 becaouse we want to start numerating nodes from 0
		nodes.push_back(mynode);
		if (y >= H)
		{
			y = 0.0;
			y -= dy;
			x += dx;
		}
		y += dy;
		if (i <= nN - nL)
		{
			if (i%nL != 0)
			{
				std::array<unsigned int, 4> IDs;
				IDs[0] = i;
				IDs[1] = i + nH;
				IDs[2] = i + nH + 1;
				IDs[3] = i + 1;
				Element myNewElement(IDs);
				elements.push_back(myNewElement);
			}
		}
	}
}
void Grid::setGlobalMatrixH()
{
	const size_t length = nodes.size();
	unsigned int row, column;
	DbMtrx H;

	this->GlobalH.resize(length);
	for (auto &x : GlobalH)
	{
		x.resize(length);
	}
	for (auto &x : elements)
	{
		H = x.getMatrixH();
		for (size_t i = 0; i < 4; i++)
		{
			for (size_t	j = 0; j < 4; j++)
			{
				row = x.getID(i)-1;
				column = x.getID(j)-1;
				GlobalH[row][column] += H[i][j];
			}
		}
	}
}
void Grid::setGlobalMatrixC()
{
	const size_t length = nodes.size();
	unsigned int row, column;
	DbMtrx C;

	this->GlobalC.resize(length);
	for (auto &x : GlobalC)
	{
		x.resize(length);
	}
	for (auto &x : elements)
	{
		C = x.getMatrixC();
		for (size_t i = 0; i < 4; i++)
		{
			for (size_t j = 0; j < 4; j++)
			{
				row = x.getID(i) - 1;
				column = x.getID(j) - 1;
				GlobalC[row][column] += C[i][j];
			}
		}
	}
	for (size_t i = 0; i < GlobalC.size(); i++)
	{
		for (size_t j = 0; j < GlobalC[0].size(); j++)
		{
			GlobalC[i][j] /= Tau;
		}
	}
}
void Grid::setGlobalMatrixHC()
{
	for (size_t i = 0; i < GlobalH.size(); i++)
	{
		for (size_t j = 0; j < GlobalH[0].size(); j++)
		{
			GlobalH[i][j] += GlobalC[i][j];
		}
	}
}
void Grid::setGlobalP()
{
	std::vector<double> sums;
	sums.resize(nodes.size());
	GlobalP.resize(nodes.size());
	for (auto &z : nodes)
	{
		std::cout << z.gett() << "   ";
	}
	std::cout << "\n";
	for (size_t i = 0; i < GlobalC.size(); i++)
	{
		for (size_t j = 0; j < GlobalC[0].size(); j++)
		{
			sums[i] += GlobalC[i][j] * nodes[j].gett();
		}
	}
	int j;
	std::fill(GlobalP.begin(), GlobalP.end(), 0);
	if(Tau==50)
	{
		std::vector<double> localP;
		for (auto &x : elements)
		{
			localP = x.getvectorP();
			for (size_t i = 0; i < 4; i++)
			{
				j = x.getID(i) - 1;
				GlobalP[j] += localP[i];
			}
		}
	}	
	for (size_t i = 0; i < GlobalP.size(); i++)
	{
		GlobalP[i] += sums[i];
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
void Grid::CalculateTemp()
{
	
	for (size_t i = 0; i < 10; i++)
	{
		setGlobalP();
		const size_t length = 16;

		double **matrix = new double*[length];
		for (size_t i = 0; i < length; i++)
		{
			matrix[i] = new double[length + 1];
		}
		for (size_t i = 0; i < length; i++)
		{
			for (size_t j = 0; j < length + 1; j++)
			{
				if (j != length)
				{
					matrix[i][j] = GlobalH[i][j];
				}
				else
				{
					matrix[i][j] = GlobalP[i];
				}
			}
		}
		JacobiMethod(matrix, length);
	}
}
void Grid::CalculateH()
{
	for (auto &x : elements)
	{
		std::array<bool, 4> BC{};
		for (int k = 0; k < 4; k++)
		{
			int j = x.getID(k) - 1;
			if (nodes[j].getX() == 0 || nodes[j].getX() >= L)
			{
				BC[k] = 1;
			}
			else if (nodes[j].getY() == 0 || nodes[j].getY() >= H)
			{
				BC[k] = 1;
			}
		}
		x.setBC(BC);
	}

	for (auto &x : elements)
	{
		std::array<double, 4> cordsx;
		std::array<double, 4> cordsy;

		for (size_t i = 0; i < 4; ++i) // we need to get cords (x,y) for our element
		{
			unsigned int tmp = x.getID(i) - 1; // -1 becaouse ID-s are 1-4 no 0-3 :(
			cordsx[i] = nodes[tmp].getX();
			cordsy[i] = nodes[tmp].getY();
			
		}
		x.calculH(cordsx, cordsy);
	}
}
void Grid::JacobiMethod(double** matrix,const size_t length)
{
	double *N = new double[length];
	double *x = new double[length](); 
	for (size_t i = 0; i < length; i++)
	{
		N[i] = pow(matrix[i][i], -1);
		matrix[i][i] = 0; 
		for (size_t j = 0; j < length; j++)
		{
			matrix[i][j] = -1.0*(N[i] * matrix[i][j]);
		}
	}
	for (size_t k = 0; k < 15; k++)
	{
		double *tmpX = new double[length]();
		for (size_t i = 0; i < length; i++)
		{
			double Mx = 0;
			for (size_t j = 0; j < length; j++)
			{
				Mx += matrix[i][j] * x[j];
			}
			tmpX[i] = Mx + N[i] * matrix[i][length];
			Mx = 0;
		}
		for (size_t i = 0; i < length; i++)
		{
			x[i] = tmpX[i];
		}
	}

	for (size_t i = 0; i < 16; i++)
	{
		std::cout << "x" << i << ": " << x[i] << "\n";
		nodes[i].sett(x[i]);
	}
	std::cout << "\n";
}
