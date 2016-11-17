#pragma once
#include "Element.h"
#include "LinearSystem.h"
#include <vector>
#include <fstream>

class Mesh
{
	void logToFile(std::fstream&, std::vector<double>&);
public:
	double q, alpha, ambientTemp;
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