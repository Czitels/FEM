#pragma once
#include "node.h"
#include <vector>
#include <array>
class Element
{
public:
	Element(std::array<unsigned int, 4>);
	void calculH(std::array<double,4>&,std::array<double,4>&);
	unsigned int getID(unsigned int);
	std::array<unsigned int, 4>& getArray();

private:
	using DbMtrx = std::vector<std::vector<double>>;
	std::array<unsigned int, 4> ID;
	DbMtrx H;
};