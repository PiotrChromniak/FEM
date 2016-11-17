#pragma once
#include "Node.h"
#include <math.h>
#include <array>

class Mesh;

class Element
{	
private:
	Mesh& _mesh;
public:
	Node &node_0, &node_1;
	const double k, S;
	std::array<std::array<double, 2>, 2> H;
	std::array<double, 2> P;

	double L() const;
	void calculateLocalMatrix();
	Element(Mesh&, Node&, Node&, double, double);
};