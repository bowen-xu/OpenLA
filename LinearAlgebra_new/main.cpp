#include <iostream>
#include "main.h"
#include <stdio.h>

using namespace std;

int main()
{
	Matrix<float> mat1(3, 3),  mat;
	
	mat(3, 3);
	mat =
		1, 4, 9,
		2, 7, 7,
		9, 6, 5;

	mat.print();

	printf("\n––¡– Ω£∫%.3f\n\n", mat.Determinant());


	return 1;
}