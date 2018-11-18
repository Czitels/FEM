#pragma once
#include "element.h"
#include "node.h"
#include <vector>
class Grid
{
public:
	Grid();
	~Grid();
	void setGrid();
	void CalculateH();
	void printGrid();

private:
	using DbMtrx = std::vector<std::vector<double>>;
	float H;
	float L;
	size_t nH;
	size_t nL;
	float t0;
	size_t nN;
	size_t nE;
	std::vector<Node> nodes;
	std::vector<Element> elements;
	DbMtrx DerivKsi;
	DbMtrx DerivEta;
};
