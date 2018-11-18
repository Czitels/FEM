#include "universal_element.h"
#include <iostream>

Universal_element::Universal_element()
{
	DerivEta.resize(4);
	DerivKsi.resize(4);
	for (size_t i = 0; i < 4; i++)
	{
		std::vector<double> row;
		row.push_back(0.25*(1 - ksi[i])*(1 - eta[i]));
		row.push_back(0.25*(1 + ksi[i])*(1 - eta[i]));
		row.push_back(0.25*(1 + ksi[i])*(1 + eta[i]));
		row.push_back(0.25*(1 - ksi[i])*(1 + eta[i]));
		Ni.push_back(row);
	}
	for (size_t i = 0; i < eta.size(); i++)
	{
		DerivEta[0].push_back(-0.25*(1 - eta[i]));
		DerivEta[1].push_back(0.25*(1 - eta[i]));
		DerivEta[2].push_back(0.25*(1 + eta[i]));
		DerivEta[3].push_back(-0.25*(1 + eta[i]));
	}
	for (size_t i = 0; i < ksi.size(); i++)
	{
		DerivKsi[0].push_back(-0.25*(1 - ksi[i]));
		DerivKsi[1].push_back(-0.25*(1 + ksi[i]));
		DerivKsi[2].push_back(0.25*(1 + ksi[i]));
		DerivKsi[3].push_back(0.25*(1 - ksi[i]));
	}
}
