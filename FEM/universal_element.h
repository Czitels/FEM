#pragma once
#include <array>
#include <vector>
struct Universal_element
{
	Universal_element();
	using DbMtrx = std::vector<std::vector<double>>;

	const double ro = 7800;
	const double c = 700;
	const double k = 25;
	const double alfa = 300;
	const double tinf = 1200;
	DbMtrx Ni; // shape functions
	DbMtrx DerivKsi;
	DbMtrx DerivEta;
};