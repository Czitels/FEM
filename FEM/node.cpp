#include "node.h"
#include <iostream>

float Node::getX()
{
	return this->x;
}
float Node::getY()
{
	return this->y;
}
float Node::gett()
{
	return this->t;
}
void Node::setBC(bool BC)
{
	this->BC = BC;
}
void Node::sett(float t)
{
	this->t = t;
}
bool Node::getBC()
{
	return this->BC;
}