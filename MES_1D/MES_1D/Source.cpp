#include "Mesh.h"
#include <chrono>

using namespace std;

int main() { 
	try {
		Mesh n;
		auto start = chrono::steady_clock::now();
		n.loadFromFile("data.txt");
		n.calculateLocalMatrices();
		n.mergeMatrices();
		n.proceedSolving();
		auto stop = chrono::steady_clock::now();
		cout << "Process took " << chrono::duration_cast<chrono::microseconds>(stop - start).count() << " microseconds\n";
	}
	catch (exception e) {
		cerr << e.what() << "\n";
	}
	system("pause");
	return 0;
}