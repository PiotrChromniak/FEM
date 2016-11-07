#include "Element.h"
#include "Mesh.h"

double Element::L() const
{ return fabs(node_1.x - node_0.x); }

Element::Element(Mesh& ref_mesh, Node &ref_1, Node &ref_2, double K) : _mesh(ref_mesh), node_0(ref_1), node_1(ref_2), k(K), H{ 0,0,0,0 }, P{0,0}
{}

void Element::calculateLocalMatrix()
{
	double temp_L = L();
	H[0][0] = (_mesh.S * k) / temp_L;
	H[1][1] = H[0][0];
	H[0][1] = -(_mesh.S * k) / temp_L;
	H[1][0] = H[0][1];

	if (NodeType::INTER != node_0.getNodeType()) {
		if (NodeType::IN == node_0.getNodeType()) P[0] = _mesh.q * _mesh.S;
		else {
			P[0] = -_mesh.alpha *  _mesh.ambientTemp *  _mesh.S;
			H[0][0] += _mesh.alpha *  _mesh.S;
		}
	}
	
	if (NodeType::INTER != node_1.getNodeType()) {
		if (NodeType::IN == node_1.getNodeType()) P[1] = _mesh.q *  _mesh.S;
		else {
			P[1] = -_mesh.alpha *  _mesh.ambientTemp *  _mesh.S;
			H[1][1] += _mesh.alpha *  _mesh.S;
		}
	}
}