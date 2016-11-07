#include "Node.h"

NodeType Node::getNodeType() const
{ return type_; }

Node::Node(int id, double X, NodeType nType) : x(X), type_(nType), ID(id)
{}