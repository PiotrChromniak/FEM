#include "Element.h"
#include "Mesh.h"

double Element::L() const
{ return fabs(node_1.x - node_0.x); }

Element::Element(Mesh& ref_mesh, Node &ref_1, Node &ref_2, double K) : mesh_(ref_mesh), node_0(ref_1), node_1(ref_2), k(K), H{ 0,0,0,0 }, P{0,0}
{}

void Element::calculateLocalMatrix()
{
	double temp_L = L();
	H[0][0] = (mesh_.S * k) / temp_L;
	H[1][1] = H[0][0];
	H[0][1] = -(mesh_.S * k) / temp_L;
	H[1][0] = H[0][1];

	if (NodeType::INTER != node_0.getNodeType()) {
		if (NodeType::IN == node_0.getNodeType()) P[0] = mesh_.q * mesh_.S;
		else {
			P[0] = -mesh_.alpha *  mesh_.ambientTemp *  mesh_.S;
			H[0][0] += mesh_.alpha *  mesh_.S;
		}
	}
	
	if (NodeType::INTER != node_1.getNodeType()) {
		if (NodeType::IN == node_1.getNodeType()) P[1] = mesh_.q *  mesh_.S;
		else {
			P[1] = -mesh_.alpha *  mesh_.ambientTemp *  mesh_.S;
			H[1][1] += mesh_.alpha *  mesh_.S;
		}
	}
}