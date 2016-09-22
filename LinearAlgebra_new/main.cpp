#include <iostream>
#include "main.h"
#include <stdio.h>

using namespace std;

int main()
{
	Matrix<float> Coff(2, 2), constant(2, 1), solution(2, 1);
	//�ⷽ����
	// 2*x +   y = 1
	// x   + 5*y = 0
	//����� x = 5/9 = 0.56, y = -1/9 = -0.11;

	Coff =			// ϵ��
		2, 1,
		1, 5;
	constant =		// ������
		1,
		0;

	//����� ���� ϵ��������� ��� �������� 
	solution = Coff.Inverse() * constant;

	solution.print();

	return 1;
}
