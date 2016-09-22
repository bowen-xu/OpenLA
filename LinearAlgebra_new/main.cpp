#include <iostream>
#include "main.h"
#include <stdio.h>

using namespace std;

int main()
{
	Matrix<float> Coff(2, 2), constant(2, 1), solution(2, 1);
	//解方程组
	// 2*x +   y = 1
	// x   + 5*y = 0
	//结果： x = 5/9 = 0.56, y = -1/9 = -0.11;

	Coff =			// 系数
		2, 1,
		1, 5;
	constant =		// 常数项
		1,
		0;

	//解矩阵 等于 系数矩阵的逆 左乘 常数矩阵 
	solution = Coff.Inverse() * constant;

	solution.print();

	return 1;
}
