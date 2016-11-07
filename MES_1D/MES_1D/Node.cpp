#include "Node.h"

NodeType Node::getNodeType() const
{ return _type; }

Node::Node(int id, double X, NodeType nType) : x(X), _type(nType), ID(id)
{}