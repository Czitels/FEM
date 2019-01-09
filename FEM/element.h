#pragma once
#include "node.h"
#include "universal_element.h"
#include <vector>
#include <array>
class Element
{
	using DbMtrx = std::vector<std::vector<double>>;

public:
	Element(std::array<unsigned int, 4>);
	void calculH(std::array<double,4>&,std::array<double,4>&);
	void calculC(std::vector<double>&, Universal_element&);
	void calculHBC(std::array<double, 4>&, std::array<double, 4>&, Universal_element&);
	unsigned int getID(unsigned int);
	void setBC(std::array<bool, 4> BC);
	DbMtrx &getMatrixH();
	DbMtrx &getMatrixC();
	std::vector<double>& getvectorP();
	std::array<unsigned int, 4>& getArray();

private:
	std::array<bool, 4> BC;
	std::array<unsigned int, 4> ID;
	std::vector<double> P;
	DbMtrx H;
	DbMtrx C;
	DbMtrx bcH;
};