#pragma once
//------------------------------------------------------------------------
class ControlOMP
{
	// Заводится один глобальный объект этого класса на всю программу и
	// управление числом потоков в OpenMP осуществляется через него
private:
	bool isInit;
	int NUMBER_OF_THREADS; // текущее число потоков в OpenMP
	int MAX_THREADS; // максимально возможное число потоков в OpenMP
	int nMinDotProduct;
	int nMinSparseMultMV;
	int nMinSparseMultMV2x2;
public:

	void InitNumThreads(); // вызывается один раз в самом начале программы
	void SetNumThreads(int num);
	int GetNumberOfThreads();
	int GetMaxThreads();
	int GetNMinDotProduct();
	int GetNMinSparseMultMV();
	int GetNMinSparseMultMV2x2();


	ControlOMP();
	~ControlOMP();
};

