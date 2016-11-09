#include "Mesh.h"


void Mesh::loadFromFile(const char *fileName)
{
	std::fstream file;
	file.open(fileName, std::ios::in);
	if (file.good()) {
		int nNode;
		file >> nNode;
		nodeVect.reserve(nNode);

		double t_x;
		int t_ID0, t_ID1, t_type;

		/* Wczytywanie node'�w */
		for (int i = 0; i < nNode; ++i) {
			file >> t_ID0 >> t_x >> t_type;
			NodeType t_nodeType = [&]() {
				switch (t_type) {
					case 0:
						return NodeType::IN;
					case 2:
						return NodeType::OUT;
					default:
						return NodeType::INTER;
				}
			}();

			nodeVect.push_back(Node(t_ID0, t_x, t_nodeType));
		}

		/* Wczytywanie element�w - korzystaj� z referencji do juz istniej�cych node'�w w nodeVec[] */
		file >> nNode;						// tak naprawd� to nElement
		elemVect.reserve(nNode);
		for (int i = 0; i < nNode; ++i) {
			file >> t_ID0 >> t_ID1 >> t_x;		// t_x to tutaj wsp�czynnik k elementu
			elemVect.push_back(Element(*this, nodeVect[t_ID0], nodeVect[t_ID1], t_x));
		}

		/* Wczytywanie warunk�w brzegowych */
		file >> q >> S >> alpha >> ambientTemp;

		file.close();
	}
	else throw 1;
}

void Mesh::calculateLocalMatrices()
{
	for (auto& element : elemVect) element.calculateLocalMatrix();
}

void Mesh::mergeMatrices()
{
	size_t size = elemVect.size();
	H.resize(size + 1);
	P.resize(size + 1);

	for (auto& H_row : H) H_row.resize(size + 2);

	for (const auto& elem : elemVect) {
		H[elem.node_0.ID][elem.node_0.ID] += elem.H[0][0];
		H[elem.node_0.ID][elem.node_1.ID] += elem.H[0][1];
		H[elem.node_1.ID][elem.node_0.ID] += elem.H[1][0];
		H[elem.node_1.ID][elem.node_1.ID] += elem.H[1][1];

		P[elem.node_0.ID] += elem.P[0];
		P[elem.node_1.ID] += elem.P[1];
	}


	/* Przygotowanie macierzy do oblicze� w bibliotece LinearSystem */
	//std::vector<double>::iterator P_it = P.begin(), migrate_it;
	//for (auto H_it = H.begin(); P_it != P.end();  ++H_it, ++P_it) {
	//	migrate_it = H_it->end() - 1;
	//	*migrate_it = -(*P_it);
	//}
	auto P_it = P.cbegin();
	std::vector<double>::iterator lastInRow;
	for (auto& row : H) {
		lastInRow = row.end() - 1;
		*lastInRow = -(*P_it);
		++P_it; 
	}
}

void Mesh::proceedSolving()
{
	/* zamienianie macierzy i wektora obci��e� z obiektem LinearSystem solver */
	//solver.swapResources(H, P);
	//solver.solveGauss();
	//std::vector<double> ans = solver.getAnsVector();
	//solver.swapResources(H,P);

	auto ans = LinearSystem::solveGauss(H, P);
	std::cout << "Metoda eliminacji Gaussa:\n";
	for (const auto& node : nodeVect)
		std::cout << "\nt" << node.ID << " = " << ans[node.ID];
	std::cout << '\n';
}

void Mesh::proceedSolvingIterational(int iterations)
{
	//solver.swapResources(H, P);
	//solver.setZeroStartVec();
	//solver.solveJacobi(iterations);
	//std::vector<double> ans = solver.getAnsVector();
	//solver.swapResources(H, P);
}