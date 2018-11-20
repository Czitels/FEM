#include <iostream>
#include <vector>
#include "universal_element.h"
#include "grid.h"

std::ostream& operator<<(std::ostream&, const Node&);
int main()
{
	Grid mygrid;
	mygrid.setGrid();
	mygrid.printGrid();
	//mygrid.CalculateH();

	system("pause");
	return 0;
}
std::ostream& operator<<(std::ostream& os, const Node& mynode)
{
	os << "ID: " << mynode.ID << " x: " << mynode.x << " y: " << mynode.y << " T: " << mynode.t;
	return os;
}
