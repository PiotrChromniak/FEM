#include "LinearSystem.h"

void LinearSystem::loadGaussMatrix()
{
	std::fstream file;
	file.open("gauss.txt", std::ios::in);
	if (file.good()) {
		size_t N;
		double a;
		file >> N;
		
		if (M.size() > 0)	M.clear();
		if (b.size() > 0)	b.clear();
		M.reserve(N);
		for (size_t i = 0; i < N; ++i)
			M.push_back(std::vector<double>(N + 1));	//N+1 za wczasu (i tak rozszerzamy macierz o wektor X[])

		//Wype³nianie wektora
		b.reserve(N);
		for (size_t i = 0; i < N; ++i) {
			file >> a;
			b.push_back(a);
		}

		//Wype³nianie macierzy
		for (size_t i = 0; i < N; ++i) {
			for (unsigned int j = 0; j < N; ++j) {
				file >> a;
				M[i][j] = a;
			}
		}

		//Rozszerzanie macierzy o wektor X
		for (size_t i = 0; i < N; ++i)
			M[i][N] = b[i];
	}
	else throw std::runtime_error("Problem while loading Gauss Matrix");
}

void LinearSystem::loadJacobiMatrix()
{
	std::fstream file;
	file.open("jacobi.txt", std::ios::in);
	if (file.good()) {
		unsigned int N;
		double a;
		file >> N;

		if (M.size() > 0)	M.clear();
		if (b.size() > 0)	b.clear();
		if (Xn.size() > 0)	Xn.clear();
		M.reserve(N);
		for (unsigned int i = 0; i < N; ++i)
			M.push_back(std::vector<double>(N));	

		//Wypelnianie wektora b
		b.reserve(N);
		for (unsigned int i = 0; i < N; ++i) {
			file >> a;
			b.push_back(a);
		}

		//Wype³nianie wekt. warunku poczatkowego
		Xn.reserve(N);
		for (unsigned int i = 0; i < N; ++i) {
			file >> a;
			Xn.push_back(a);
		}

		//Wype³nianie macierzy
		for (unsigned int i = 0; i < N; ++i) {
			for (unsigned int j = 0; j < N; ++j) {
				file >> a;
				M[i][j] = -a;
			}
		}
	}
	else throw std::runtime_error("Problem while loading Jacobi Matrix");
}

void LinearSystem::solveGauss()
{
	/* Zerowanie macierzy */
	for (unsigned int k = 0; k < M.size() + 1; k++) {
		for (unsigned int i = 1 + k; i < M.size(); i++)
		{
			if (fabs(M[k][k] < eps_)) throw 1.0;
			double coefficient = M[i][k] / M[k][k];
			for (unsigned int j = k; j < M[i].size(); j++)
				M[i][j] -= coefficient * M[k][j];
		}
	}

	/* Redukowanie macierzy o wektor X */
	int last = b.size();
	for (unsigned int i = 0; i < b.size(); ++i) {
		b[i] = M[i][last];
		M[i].pop_back();
	}

	//Wyliczanie xn-1....x0
	std::vector<double> X(b.size());
	int n = b.size() - 1;
	//if (n == b.size() - 1)
	X[n] = b[n] / M[n][n];
	if (b.size() > 1) {
		int k = 2;
		for (int i = b.size() - 2; i >= 0; i--, k++) {

			double eq = b[i];
			for (int j = b.size() - 1; j > b.size() - k; j--) {
				eq -= X[j] * M[i][j];
			}
			if (fabs(M[i][i] < eps_)) throw 1.0;
			X[i] = eq / M[i][i];
		}
	}

	//std::cout << "Metoda eliminacji Gaussa:\n";
	//for (unsigned int i = 0; i < X.size(); ++i)
	//	std::cout << "\nx" << i << "= " << X[i];
	//std::cout <<std:: endl;
	Xn.swap(X);
}

void LinearSystem::solveJacobi(int steps)
{
	checkDiagonalDominance();
	/* invD to odwrotnoœæ macierzy D */
	std::vector<std::vector<double>> invD(M.size(), std::vector<double>(M.size()));
	for (unsigned int i = 0; i < M.size(); ++i) {
		invD[i][i] = -1.0 / M[i][i];
		M[i][i] = 0.0;
	}

	//  M_0 = (U+L) x invD - oszczêdzamy sobie mno¿enia na pizdeczke za ka¿dym razem 
	//	V_0 = invD x b jak wy¿ej 
	std::vector<std::vector<double>> M_0(M.size(), std::vector<double>(M.size()));
	std::vector<double> V_0(M.size());
	unsigned int N = M.size(); 

	/* Korzystam z faktu ze macierz diagonalna ma wartosci tylko na przekatnej - oszczednosc mnozenia/dodawania = profit :^) */
	for (unsigned int y = 0; y < N; ++y) 
		for (unsigned int x = 0; x < N; ++x) 
			M_0[y][x] = invD[y][y] * M[y][x];

	//Macierz M ju¿ nie jest potrzebna
	M.clear();

	/* Znowu korzystamy z faktu macierzy diagonalnej - znowu profit :^) */
	for (unsigned int y = 0; y < N; ++y)
		V_0[y] = b[y] * invD[y][y];
	b.clear();

	/* Iteracyjne generowanie rozwiazan - im wiecej iteracji tym bardziej zblizony wynik (patrz wynik z Gaussa) */
	std::cout << "\nMetoda Jacobiego ("<<steps<<" krokow):\n\nx0" << " = [";
	for (unsigned int x = 0; x < N; ++x) {
		std::cout << Xn[x];
		if (x != N - 1) std::cout << ", ";
	}
	std::cout << "]"; 

	std::vector<double> temp(N);
	for (int i = 1; i < steps + 1; ++i) {
		for (unsigned int y = 0; y < N; ++y) {
			double sum = 0.0;
			for (unsigned int j = 0; j < N; ++j)
				sum += Xn[j] * M_0[y][j];
			sum += V_0[y];
			temp[y] = sum;
		}
		std::cout << "\nx" << i << " = [";
		for (unsigned int x = 0; x < N; ++x) {
			std::cout << temp[x];
			if (x != N - 1) std::cout << ", ";
		}
		std::cout << "]";
		Xn = temp;
	}
	std::cout << '\n';
}

void LinearSystem::setZeroStartVec()
{
	Xn.resize(b.size());
}

std::vector<double> LinearSystem::tmpJacobi(std::vector<std::vector<double>> &M, std::vector<double> &b, std::vector<double>& Xn, int steps)
{
	checkDiagonal(M);
	/* invD to odwrotnoœæ macierzy D */
	std::vector<std::vector<double>> invD(M.size(), std::vector<double>(M.size()));
	for (unsigned int i = 0; i < M.size(); ++i) {
		invD[i][i] = -1.0 / M[i][i];
		M[i][i] = 0.0;
	}

	//  M_0 = (U+L) x invD - oszczêdzamy sobie mno¿enia na pizdeczke za ka¿dym razem 
	//	V_0 = invD x b jak wy¿ej 
	std::vector<std::vector<double>> M_0(M.size(), std::vector<double>(M.size()));
	std::vector<double> V_0(M.size());
	unsigned int N = M.size();

	/* Korzystam z faktu ze macierz diagonalna ma wartosci tylko na przekatnej - oszczednosc mnozenia/dodawania = profit :^) */
	for (unsigned int y = 0; y < N; ++y)
		for (unsigned int x = 0; x < N; ++x)
			M_0[y][x] = invD[y][y] * M[y][x];

	//Macierz M ju¿ nie jest potrzebna
	M.clear();

	/* Znowu korzystamy z faktu macierzy diagonalnej - znowu profit :^) */
	for (unsigned int y = 0; y < N; ++y)
		V_0[y] = b[y] * invD[y][y];
	b.clear();

	/* Iteracyjne generowanie rozwiazan - im wiecej iteracji tym bardziej zblizony wynik (patrz wynik z Gaussa) */
	std::cout << "\nMetoda Jacobiego (" << steps << " krokow):\n\nx0" << " = [";
	for (unsigned int x = 0; x < N; ++x) {
		std::cout << Xn[x];
		if (x != N - 1) std::cout << ", ";
	}
	std::cout << "]";

	std::vector<double> temp(N);
	for (int i = 1; i < steps + 1; ++i) {
		for (unsigned int y = 0; y < N; ++y) {
			double sum = 0.0;
			for (unsigned int j = 0; j < N; ++j)
				sum += Xn[j] * M_0[y][j];
			sum += V_0[y];
			temp[y] = sum;
		}
		std::cout << "\nx" << i << " = [";
		for (unsigned int x = 0; x < N; ++x) {
			std::cout << temp[x];
			if (x != N - 1) std::cout << ", ";
		}
		std::cout << "]";
		Xn = temp;
	}
	std::cout << '\n';

	return Xn;
}

std::vector<double> LinearSystem::tmpGauss(std::vector<std::vector<double>> &M, std::vector<double> &b, double eps)
{
	/* Zerowanie macierzy */
	for (unsigned int k = 0; k < M.size() + 1; k++) {
		for (unsigned int i = 1 + k; i < M.size(); i++)
		{
			if (fabs(M[k][k] < eps)) throw std::logic_error("Division by zero in Gauss Elimination");
			double coefficient = M[i][k] / M[k][k];
			for (unsigned int j = k; j < M[i].size(); j++)
				M[i][j] -= coefficient * M[k][j];
		}
	}

	/* Redukowanie macierzy o wektor X */
	int last = b.size();
	for (unsigned int i = 0; i < b.size(); ++i) {
		b[i] = M[i][last];
		M[i].pop_back();
	}

	//Wyliczanie xn-1....x0
	std::vector<double> X(b.size());
	int n = b.size() - 1;
	//if (n == b.size() - 1)
	X[n] = b[n] / M[n][n];
	if (b.size() > 1) {
		int k = 2;
		for (int i = b.size() - 2; i >= 0; i--, k++) {

			double eq = b[i];
			for (int j = b.size() - 1; j > b.size() - k; j--) {
				eq -= X[j] * M[i][j];
			}
			if (fabs(M[i][i] < eps)) throw std::logic_error("Division by zero in Gauss Elimination");
			X[i] = eq / M[i][i];
		}
	}

	//std::cout << "Metoda eliminacji Gaussa:\n";
	//for (unsigned int i = 0; i < X.size(); ++i)
	//	std::cout << "\nx" << i << "= " << X[i];
	//std::cout <<std:: endl;

	return X;
}

std::vector<double> LinearSystem::solveJacobiFromFile(const char *, int)
{
	return std::vector<double>();
}

void LinearSystem::checkDiagonal(const std::vector<std::vector<double>> &M)
{
	bool strictInequalityConditon = false;
	for (unsigned int i = 0; i < M.size(); ++i) {
		double rowSum = 0.0;
		for (unsigned int j = 0; j < M[i].size(); ++j) {
			if (i != j)	rowSum += fabs(M[i][j]);
		}
		if (rowSum > fabs(M[i][i]))	throw std::runtime_error("Matrix isn't diagonaly dominant");
		if (fabs(M[i][i]) > rowSum && !strictInequalityConditon)	strictInequalityConditon = true;
	}
	if (!strictInequalityConditon)	throw std::runtime_error("Matrix isn't diagonaly dominant");
	std::cout << "\nMacierz M jest diagonalnie dominujaca.\n";
}

void LinearSystem::checkDiagonalDominance()
{
	bool strictInequalityConditon = false;
	for (unsigned int i = 0; i < M.size(); ++i) {
		double rowSum = 0.0;
		for (unsigned int j = 0; j < M[i].size(); ++j) {
			if (i != j)	rowSum += fabs(M[i][j]);
		}
		if (rowSum > fabs(M[i][i]))	throw std::runtime_error("Matrix isn't diagonaly dominant");
		if (fabs(M[i][i]) > rowSum && !strictInequalityConditon)	strictInequalityConditon = true;
	}
	if (!strictInequalityConditon)	throw std::runtime_error("Matrix isn't diagonaly dominant");
	std::cout << "\nMacierz M jest diagonalnie dominujaca.\n";
}

LinearSystem::LinearSystem() : eps_(1e-12)
{}