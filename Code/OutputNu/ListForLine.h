#pragma once

struct NodeForLine
{
	long material;
	double val[3][2]; // ��� ���� ������������ ��� 3 ����������, ��� ����� - ������ ����
};
//------------------------------------------------------------------------
class ListOfNodesForLine
{
public:
	long kuzlov;
	long size; // ����� ����� � ������������ 
	vector< vector<NodeForLine> > item;

	ListOfNodesForLine();
	~ListOfNodesForLine();
	void Init(long kuzlov);
	void Add(long node, long mat);

	void CalcSize();
};
