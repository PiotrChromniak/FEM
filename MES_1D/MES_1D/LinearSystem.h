#pragma once

#include <vector>
#include <fstream>
#include <iostream>
#include <math.h>

class LinearSystem
{
	static void checkDiagonal(const std::vector<std::vector<double>>&);
	static void loadJacobiMatrix(std::fstream&, std::vector<std::vector<double>>&, std::vector<double>&, std::vector<double>&);
	static void loadGaussMatrix(std::fstream&, std::vector<std::vector<double>>&, std::vector<double>&);
public:
	static std::vector<double> solveJacobi(std::vector<std::vector<double>>&, std::vector<double>&, std::vector<double>&, int = 100);
	static std::vector<double> solveGauss(std::vector<std::vector<double>>&, std::vector<double>&, double = 1.0e-12);
	static std::vector<double> solveJacobiFromFile(const char*, int = 100);
	static std::vector<double> solveGaussFromFile(const char*, double = 1.0e-12);
	
};