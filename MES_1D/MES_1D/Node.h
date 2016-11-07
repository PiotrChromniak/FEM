#pragma once
#include "NodeType.h"

class Node
{
	const NodeType _type;
public:
	const double x;
	const int ID;

	NodeType getNodeType() const;

	Node(int, double, NodeType);
};