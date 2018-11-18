#pragma once
#include <iostream>
class Node
{
public:
	size_t ID;
	float x, y, t;
	friend std::ostream& operator<<(std::ostream&, const Node&);

	Node(size_t mID, float mx, float my, float mt):ID(mID),x(mx),y(my),t(mt){}
	float getX();
	float getY();
};