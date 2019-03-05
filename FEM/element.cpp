#include <vector>
#include <cmath>
#include "element.h"
#include "universal_element.h"

Element::Element(std::array<unsigned int, 4> IDs) :ID(IDs) { H.resize(4); }

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

	calculC(detJ, univ);
	calculHBC(cordsx, cordsy, univ);

	for (size_t i = 0; i < 4; i++)
	{
		JacMinusOne[0].push_back(Jacobian[3][i] / detJ[i]);
		JacMinusOne[1].push_back(-Jacobian[1][i] / detJ[i]);
		JacMinusOne[2].push_back(-Jacobian[2][i] / detJ[i]);
		JacMinusOne[3].push_back(Jacobian[0][i] / detJ[i]);
	}
	
	for (size_t i = 0; i < 4; i++)
	{
		for (size_t j = 0; j < 4; j++)
		{
			DerivX[i].push_back(JacMinusOne[0][i] * univ.DerivEta[j][i] + JacMinusOne[1][i] * univ.DerivKsi[j][i]);
			DerivY[i].push_back(JacMinusOne[2][i] * univ.DerivEta[j][i] + JacMinusOne[3][i] * univ.DerivKsi[j][i]);
		}
	}

	DbMtrx DerivXTransp1;
	DbMtrx DerivXTransp2;
	DbMtrx DerivXTransp3;
	DbMtrx DerivXTransp4;

	DbMtrx DerivYTransp1;
	DbMtrx DerivYTransp2;
	DbMtrx DerivYTransp3;
	DbMtrx DerivYTransp4;

	DerivXTransp1.resize(4);
	DerivXTransp2.resize(4);
	DerivXTransp3.resize(4);
	DerivXTransp4.resize(4);

	DerivYTransp1.resize(4);
	DerivYTransp2.resize(4);
	DerivYTransp3.resize(4);
	DerivYTransp4.resize(4);


	for (size_t i = 0; i < 4; i++)
	{
		for (size_t j = 0; j < 4; j++)
		{
			DerivXTransp1[i].push_back(DerivX[0][j] * DerivX[0][i]);
			DerivXTransp2[i].push_back(DerivX[1][j] * DerivX[1][i]);
			DerivXTransp3[i].push_back(DerivX[2][j] * DerivX[2][i]);
			DerivXTransp4[i].push_back(DerivX[3][j] * DerivX[3][i]);

			DerivYTransp1[i].push_back(DerivY[0][j] * DerivY[0][i]);
			DerivYTransp2[i].push_back(DerivY[1][j] * DerivY[1][i]);
			DerivYTransp3[i].push_back(DerivY[2][j] * DerivY[2][i]);
			DerivYTransp4[i].push_back(DerivY[3][j] * DerivY[3][i]);
		}
	}

	for (size_t i = 0; i < 4; i++)
	{
		for (size_t j = 0; j < 4; j++)
		{
			DerivXTransp1[i][j] *= detJ[0];
			DerivXTransp2[i][j] *= detJ[1];
			DerivXTransp3[i][j] *= detJ[2];
			DerivXTransp4[i][j] *= detJ[3];

			DerivYTransp1[i][j] *= detJ[0];
			DerivYTransp2[i][j] *= detJ[1];
			DerivYTransp3[i][j] *= detJ[2];
			DerivYTransp4[i][j] *= detJ[3];
		}
	}

	DbMtrx FinalP1;
	DbMtrx FinalP2;
	DbMtrx FinalP3;
	DbMtrx FinalP4;

	FinalP1.resize(4);
	FinalP2.resize(4);
	FinalP3.resize(4);
	FinalP4.resize(4);
	for (size_t i = 0; i < 4; i++)
	{
		for (size_t j = 0; j < 4; j++)
		{
			FinalP1[i].push_back(univ.k*(DerivXTransp1[i][j] + DerivYTransp1[i][j]));
			FinalP2[i].push_back(univ.k*(DerivXTransp2[i][j] + DerivYTransp2[i][j]));
			FinalP3[i].push_back(univ.k*(DerivXTransp3[i][j] + DerivYTransp3[i][j]));
			FinalP4[i].push_back(univ.k*(DerivXTransp4[i][j] + DerivYTransp4[i][j]));
		}
	}
	for (size_t i = 0; i < 4; i++)
	{
		for (size_t j = 0; j < 4; j++)
		{
			H[i].push_back(FinalP1[i][j] + FinalP2[i][j] + FinalP3[i][j] + FinalP4[i][j]+bcH[i][j]); // we sum points and add boundary conditions H
		}
	}
}
void Element::calculC(std::vector<double> &detJ, Universal_element &univ)
{
	C.resize(4);

	DbMtrx Point1;
	DbMtrx Point2;
	DbMtrx Point3;
	DbMtrx Point4;

	Point1.resize(4);
	Point2.resize(4);
	Point3.resize(4);
	Point4.resize(4);
	detJ = { 0.00027889,0.00027889,0.00027889,0.00027889 };
	for (size_t i = 0; i < 4; i++)
	{
		for (size_t	j = 0; j < 4; j++)
		{
			Point1[i].push_back(univ.Ni[j][0] * univ.Ni[i][0] * detJ[0] * univ.c * univ.ro);
			Point2[i].push_back(univ.Ni[j][1] * univ.Ni[i][1] * detJ[1] * univ.c * univ.ro);
			Point3[i].push_back(univ.Ni[j][2] * univ.Ni[i][2] * detJ[2] * univ.c * univ.ro);
			Point4[i].push_back(univ.Ni[j][3] * univ.Ni[i][3] * detJ[3] * univ.c * univ.ro);
		}
	}
	for (size_t i = 0; i < 4; i++)
	{
		for (size_t j = 0; j < 4; j++)
		{
			C[i].push_back(Point1[i][j] + Point2[i][j] + Point3[i][j] + Point4[i][j]);
		}
	}
}
void Element::calculHBC(std::array<double, 4> &cordsx, std::array<double, 4> &cordsy, Universal_element& univ)
{
	P.resize(4);
	bcH.resize(4);
	for (size_t i = 0; i < 4; i++)
	{
		bcH[i].resize(4);
	}
	const double KSIPC1 = -1 / sqrt(3);
	const double KSIPC2 = -KSIPC1;
	double detJ, length;
	std::vector<double> N1;
	std::vector<double> N2;

	N1.push_back(0.5 * (1.0 - KSIPC1));
	N1.push_back(0.5 * (1.0 - KSIPC2));
	N2.push_back(0.5 * (KSIPC1 + 1.0));
	N2.push_back(0.5 * (KSIPC2 + 1.0));

	DbMtrx sum;
	DbMtrx boundP;

	sum.resize(2);
	sum[0].resize(2);
	sum[1].resize(2);

	boundP.resize(4);
	for (auto &x:boundP)
	{
		x.resize(4);
	}

	for (size_t i = 0; i < 4; i++)
	{
		detJ = 0;
		length = 0;
		if (i != 3)
		{
			if (BC[i] == true && BC[i+1] == true)
			{
				double tmpx = cordsx[i + 1] - cordsx[i];
				double tmpy = cordsy[i + 1] - cordsy[i];
				length = sqrt(pow(tmpx,2) + pow(tmpy,2));
				detJ = length / 2;

				double pc1;
				double pc2;
				for (size_t	j = 0; j < 2; j++)
				{
					for (size_t	k = 0; k < 2; k++)
					{
						pc1 = N1[j] * N1[k] * univ.alfa;
						pc2 = N2[j] * N2[k] * univ.alfa;
						sum[j][k] = (pc1 + pc2)*detJ;
						bcH[j + i][k + i] += sum[j][k];
					}
				}
				P[i] += univ.tinf*(sum[0][0] + sum[1][0]);
				P[i+1] += univ.tinf*(sum[0][1] + sum[1][1]);
			}
		}
		else if(BC[i] == true && BC[0] == true)
		{
			double tmpx = cordsx[i] - cordsx[0];
			double tmpy = cordsy[i] - cordsy[0];
			length = sqrt(pow(tmpx, 2) + pow(tmpy, 2));
			detJ = length / 2;

			double pc1;
			double pc2;
			for (size_t j = 0; j < 2; j++)
			{
				for (size_t k = 0; k < 2; k++)
				{
					pc1 = N1[j] * N1[k] * univ.alfa;
					pc2 = N2[j] * N2[k] * univ.alfa;
					sum[j][k] = (pc1 + pc2)*detJ;
					bcH[j * i][k * i] += sum[j][k]; // multiply, beacaouse we want to set corners of matrix
				}
			}
			P[i] += univ.tinf*(sum[0][0] + sum[1][0]);
			P[0] += univ.tinf*(sum[0][1] + sum[1][1]);
		}
	}
}
unsigned int Element::getID(unsigned int id)
{
	return this->ID[id];
}
void Element::setBC(std::array<bool, 4> BC)
{
	this->BC = BC;
}
std::vector<std::vector<double>>& Element::getMatrixH()
{
	return H;
}
std::vector<std::vector<double>>& Element::getMatrixC()
{
	return C;
}
std::vector<double>& Element::getvectorP()
{
	return P;
}
std::array<unsigned int, 4>& Element::getArray()
{
	return this->ID;
}
