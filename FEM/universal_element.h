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
//TODO: stworzyc macierz uniwersalna i podstawiæ wzory, które oblicz¹ elementy.
//Mo¿e klasa element uniwersalny bêdzie klas¹ abstrakcyjn¹ z której bêd¹ dziedziczyæ inne?
//Macierz uniwersalna ma byæ klas¹ czy ki uj xD a mo¿e to vector vectorów?