#pragma once
class Portrait_profil
{
public:
	long n;       // ����������� ����
	long n_elem;  // ����� ��������� � �����

	long *ig; // ������ ������� ����� �����
	bool mem;

	Portrait_profil(long n_elem);
	~Portrait_profil();

	void Gen_ig_for_1d_harm_problem();
};
