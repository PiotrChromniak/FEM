#pragma once
#include "NodeType.h"

class Node
{
	const NodeType type_;
public:
	const double x;
	const int ID;

	NodeType getNodeType() const;

	Node(int, double, NodeType);
};