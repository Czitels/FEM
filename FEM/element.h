#pragma once
#include "node.h"
#include <vector>
#include <array>
class Element
{
public:
	void calculH(std::array<double,4>&,std::array<double,4>&);
	unsigned int getID(unsigned int);

private:
	using DbMtrx = std::vector<std::vector<double>>;
	std::array<unsigned int, 4> ID;
	DbMtrx H;
};