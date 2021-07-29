#include "stdafx.h"
#include "ListForLine.h"
//------------------------------------------------------------------------
ListOfNodesForLine::ListOfNodesForLine()
{
	kuzlov = 0;
	size = 0;
}
//------------------------------------------------------------------------
void ListOfNodesForLine::Init(long kuzlov)
{
	this->kuzlov = kuzlov;
	item.resize(kuzlov);
}
//------------------------------------------------------------------------
ListOfNodesForLine::~ListOfNodesForLine()
{
	long i, j, sz;
	for (i=0; i<kuzlov; i++)
		if((sz = (long)item[i].size())!=0)
			for(j=0; j<sz; j++)
				item[i].clear();				
}
//------------------------------------------------------------------------
void ListOfNodesForLine::Add(long node, long mat)
{	
	long i;
	bool flag=true;
	NodeForLine t;

	for (i=0; i<(long)item[node].size(); i++)
		if(item[node][i].material==mat)
		{
			flag = false;
			break;
		}

	if (flag)
	{
		t.material = mat;
		item[node].push_back(t);
	}
}
//------------------------------------------------------------------------
void ListOfNodesForLine::CalcSize()
{
	long i;

	size = 0;
	for (i=0; i<kuzlov; i++)
		size += (long)item[i].size();
}
