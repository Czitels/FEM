#pragma once
#include "element.h"
#include "node.h"
#include <vector>
class Grid
{
	using DbMtrx = std::vector<std::vector<double>>;

public:
	Grid();
	~Grid();
	void setGrid();
	void setGlobalMatrixH();
	void setGlobalMatrixC();
	void setGlobalMatrixHC();
	void setGlobalP();
	void CalculateTemp();
	void CalculateH();
	void printGrid();
	void JacobiMethod(double** Matrix, const size_t);

private:
	float H;
	float L;
	int Tau = 50;
	size_t nH;
	size_t nL;
	float t0;
	size_t nN;
	size_t nE;
	std::vector<Node> nodes;
	std::vector<Element> elements;
	std::vector<double> GlobalP;
	DbMtrx GlobalH;
	DbMtrx GlobalC;
};
