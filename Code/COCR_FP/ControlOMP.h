#pragma once
//------------------------------------------------------------------------
class ControlOMP
{
	// ��������� ���� ���������� ������ ����� ������ �� ��� ��������� �
	// ���������� ������ ������� � OpenMP �������������� ����� ����
private:
	bool isInit;
	int NUMBER_OF_THREADS; // ������� ����� ������� � OpenMP
	int MAX_THREADS; // ����������� ��������� ����� ������� � OpenMP
	int nMinDotProduct;
	int nMinSparseMultMV;
	int nMinSparseMultMV2x2;
public:

	void InitNumThreads(); // ���������� ���� ��� � ����� ������ ���������
	void SetNumThreads(int num);
	int GetNumberOfThreads();
	int GetMaxThreads();
	int GetNMinDotProduct();
	int GetNMinSparseMultMV();
	int GetNMinSparseMultMV2x2();


	ControlOMP();
	~ControlOMP();
};

