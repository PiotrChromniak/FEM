#pragma once
#include "Element.h"
#include "LinearSystem.h"
#include <vector>
#include <fstream>

class Mesh
{
public:
	double q, S, alpha, ambientTemp;
	std::vector<Node> nodeVect;
	std::vector<Element> elemVect;
	std::vector<std::vector<double>> H;
	std::vector<double> P;

	void loadFromFile(const char*);
	void calculateLocalMatrices();
	void mergeMatrices();
	void proceedSolving();
	void proceedSolvingIterational(int = 100);
};