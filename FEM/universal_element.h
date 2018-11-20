#pragma once
#include <array>
#include <vector>
struct Universal_element 
{
	Universal_element();
	using DbMtrx = std::vector<std::vector<double>>;
	std::array<const double, 4> ksi = { -0.5774, 0.5773, 0.5773, -0.5774 };
	std::array<const double, 4> eta = { -0.5774, -0.5774, 0.5773, 0.5773 };
	DbMtrx Ni; // shape functions
	DbMtrx DerivKsi;
	DbMtrx DerivEta;
};
