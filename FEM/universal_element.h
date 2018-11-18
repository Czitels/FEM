#pragma once
#include <array>
#include <vector>
struct Universal_element //singleton?
{
	Universal_element();
	using DbMtrx = std::vector<std::vector<double>>;
	std::array<const double, 4> ksi = { -0.5774, 0.5773, 0.5773, -0.5774 };
	std::array<const double, 4> eta = { -0.5774, -0.5774, 0.5773, 0.5773 };
	DbMtrx Ni; // shape functions
	DbMtrx DerivKsi;
	DbMtrx DerivEta;
};
//TODO: stworzyc macierz uniwersalna i podstawi� wzory, kt�re oblicz� elementy.
//Mo�e klasa element uniwersalny b�dzie klas� abstrakcyjn� z kt�rej b�d� dziedziczy� inne?
//Macierz uniwersalna ma by� klas� czy ki uj xD a mo�e to vector vector�w?