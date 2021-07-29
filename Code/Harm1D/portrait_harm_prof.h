#pragma once
class Portrait_profil
{
public:
	long n;       // размерность СЛАУ
	long n_elem;  // число элементов в сетке

	long *ig; // массив адресов начал строк
	bool mem;

	Portrait_profil(long n_elem);
	~Portrait_profil();

	void Gen_ig_for_1d_harm_problem();
};
