#include "Mesh.h"


void Mesh::loadFromFile(const char *fileName)
{
	std::fstream file;
	file.open(fileName, std::ios::in);
	if (file.good()) {
		int nodeNum;
		file >> nodeNum;
		nodeVect.reserve(nodeNum);

		double t_x;
		int t_ID, t_type;

		/* Wczytywanie node'�w */
		for (int i = 0; i < nodeNum; ++i) {
			file >> t_ID >> t_x >> t_type;
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

			nodeVect.push_back(Node(t_ID, t_x, t_nodeType));
		}

		/* Wczytywanie element�w - korzystaj� z referencji do juz istniej�cych node'�w w nodeVect[] */
		int elementNum, t_ID1, t_ID2;
		double t_S;
		file >> elementNum;
		elemVect.reserve(elementNum);
		for (int i = 0; i < elementNum; ++i) {
			file >> t_ID1 >> t_ID2 >> t_x >> t_S;		// t_x to tutaj wsp�czynnik k elementu
			elemVect.push_back(Element(*this, nodeVect[t_ID1], nodeVect[t_ID2], t_x, t_S));
		}

		/* Wczytywanie warunk�w brzegowych */
		file >> q >> alpha >> ambientTemp;

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
	auto answer = LinearSystem::solveGauss(H, P);
	std::cout << "Metoda eliminacji Gaussa:\n";
	for (const auto& node : nodeVect)
		std::cout << "\nt" << node.ID << " = " << answer[node.ID];
	std::cout << '\n';
	std::fstream log("log.txt", std::ios::out);
	if (log.good())
		logToFile(log, answer);
	else
		throw std::runtime_error("Couldn't make log file.");
}

void Mesh::proceedSolvingJacobi(int iterations)
{
	/* TODO */
}

void Mesh::logToFile(std::fstream &file, std::vector<double> &ans)
{
	for (auto const& partOfVector : ans)
		file << partOfVector << '\n';
}