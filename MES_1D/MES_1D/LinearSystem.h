#pragma once

#include <vector>
#include <fstream>
#include <iostream>
#include <math.h>

class LinearSystem
{
	const double eps_;
	std::vector<std::vector<double>> M;
	std::vector<double> b,Xn;
	void checkDiagonalDominance();
public:
	LinearSystem();
	void loadGaussMatrix();
	void loadJacobiMatrix();
	void solveGauss();
	void solveJacobi(int);
	void setZeroStartVec();
	static std::vector<double> tmpJacobi(std::vector<std::vector<double>>&, std::vector<double>&, std::vector<double>&, int = 100);
	static std::vector<double> tmpGauss(std::vector<std::vector<double>>&, std::vector<double>&, double = 1.0e-12);
	static std::vector<double> solveJacobiFromFile(const char*, int = 100);
	static std::vector<double> solveGaussFromFile(const char*, double = 1.0e-12);
	static void checkDiagonal(const std::vector<std::vector<double>>&);
};