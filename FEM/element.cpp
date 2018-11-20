#include <vector>
#include "element.h"
#include "universal_element.h"

void Element::calculH(std::array<double,4> &cordsx, std::array<double,4> &cordsy)
{
	Universal_element univ;
	DbMtrx Jacobian;
	std::vector<double> detJ;
	DbMtrx JacMinusOne;
	DbMtrx DerivX;
	DbMtrx DerivY;

	Jacobian.resize(4);
	JacMinusOne.resize(4);
	DerivX.resize(4);
	DerivY.resize(4);

	for (size_t i = 0; i < 4; i++)
	{
		double valueXEta = 0, valueYEta = 0;
		double valueXKsi = 0, valueYKsi = 0;
		for (size_t j = 0; j < 4; j++)
		{
			valueXEta += univ.DerivEta[j][i] * cordsx[j];
			valueYEta += univ.DerivEta[j][i] * cordsy[j];
			valueXKsi += univ.DerivKsi[j][i] * cordsx[j];
			valueYKsi += univ.DerivKsi[j][i] * cordsy[j];
		}
		Jacobian[0].push_back(valueXEta);
		Jacobian[1].push_back(valueYEta);
		Jacobian[2].push_back(valueXKsi);
		Jacobian[3].push_back(valueYKsi);
	}

	for (size_t i = 0; i < 4; i++)
	{
		double value = Jacobian[0][i] * Jacobian[3][i] - Jacobian[1][i] * Jacobian[2][i];
		detJ.push_back(value);
	}

	for (size_t i = 0; i < 4; i++)
	{
		JacMinusOne[0].push_back(Jacobian[3][i] / detJ[i]);
		JacMinusOne[1].push_back(Jacobian[1][i] / detJ[i]);
		JacMinusOne[2].push_back(Jacobian[2][i] / detJ[i]);
		JacMinusOne[3].push_back(Jacobian[0][i] / detJ[i]);
	}

	for (size_t i = 0; i < 4; i++)
	{

	}
}

unsigned int Element::getID(unsigned int id)
{
	return this->ID[id];
}