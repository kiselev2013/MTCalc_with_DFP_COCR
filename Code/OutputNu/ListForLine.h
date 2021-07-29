#pragma once

struct NodeForLine
{
	long material;
	double val[3][2]; // для узла используются все 3 компоненты, для ребра - только одна
};
//------------------------------------------------------------------------
class ListOfNodesForLine
{
public:
	long kuzlov;
	long size; // число узлов с повторениями 
	vector< vector<NodeForLine> > item;

	ListOfNodesForLine();
	~ListOfNodesForLine();
	void Init(long kuzlov);
	void Add(long node, long mat);

	void CalcSize();
};
