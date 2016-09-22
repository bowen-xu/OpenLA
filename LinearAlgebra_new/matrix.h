#ifndef _CHAPTER1_H_
#define _CHAPTER1_H_
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

template<typename Tp>
class Matrix;
class noEqual {};
class cannotMultiply {};

enum RorCOpType_1 //行列变换类型枚举
{
	EXCHANGE_   //换法变换
};
enum RorCOpType_2 //行列变换类型枚举
{
	MULTIPLY_, //倍法变换
	ELIMINATE_, //消法变换
};


template<typename Tp>
class Matrix
{
private:
	int		rows;
	int		cols;
	Tp		*pMat;
	int		iterator;

	void	create();
	void	release();
	void	iEqual(const Matrix<Tp> &mat1, const Matrix<Tp> &mat2);
	void	iCanMultiply(const Matrix<Tp> &mat1, const Matrix<Tp> &mat2);

public:
	Matrix();
	Matrix(const int &row, const int &col);



	void		print();
	void		input();
	Matrix<Tp>	zeros();
	Matrix<Tp>	ones();
	Matrix<Tp>	identity();
	Matrix<Tp>	transpose();
	Matrix<Tp>	RowOp(RorCOpType_1 op_type, int row1, int row2);			//初等行变换
	Matrix<Tp>	RowOp(RorCOpType_2 op_type, int row1, double k);			//初等行变换
	Matrix<Tp>	RowOp(RorCOpType_2 op_type, int row1, double k, int row2);	//初等行变换
	Matrix<Tp>	ColOp(RorCOpType_1 op_type, int col1, int col2);			//初等列变换
	Matrix<Tp>	ColOp(RorCOpType_2 op_type, int col1, double k);			//初等列变换
	Matrix<Tp>	ColOp(RorCOpType_2 op_type, int col1, double k, int col2);	//初等列变换
	long double	Determinant();
	Matrix<Tp>	Inverse();

	Matrix<Tp> operator =  (Matrix<Tp> mat);
	Matrix<Tp> operator +  (const Matrix<Tp> &mat);
	Matrix<Tp> operator -  (const Matrix<Tp> &mat);
	Matrix<Tp> operator *  (const Matrix<Tp> &mat);
	Matrix<Tp> operator *  (const Tp num);
	Matrix<Tp> operator /  (const Tp num);
	void	   operator () (const int row, const int col);
	Matrix<Tp>&operator =  (const Tp ele);
	Matrix<Tp>&operator ,  (const Tp ele);
	template<typename _Tp>
	friend Matrix<_Tp> operator - (const Matrix<_Tp> &mat);
	template<typename _Tp>
	friend Matrix<_Tp> operator + (const Matrix<_Tp> &mat);
};


template<typename Tp>
Matrix<Tp>& Matrix<Tp>::operator =  (const Tp ele)
{
	iterator = 0;
	return operator ,(ele);
}

template<typename Tp>
Matrix<Tp>& Matrix<Tp>::operator ,  (const Tp ele)
{
	static const int iterator_max = rows * cols;
	if (iterator >= iterator_max)
	{
		iterator = 0;
	}
	pMat[iterator] = ele;
	iterator++;
	return *this;
}

template<typename Tp>
Matrix<Tp> Matrix<Tp>::operator = (Matrix<Tp> mat)
{
	if (rows == 0 || cols == 0 || rows != mat.rows || cols != mat.cols)
	{
		release();
		rows = mat.rows;
		cols = mat.cols;
		create();
	}
	Tp *p1 = pMat;
	Tp *p2 = mat.pMat;
	for (int i = 0; i < mat.rows; i++)
	{
		for (int j = 0; j < mat.cols; j++)
		{
			*p1 = *p2;
			p1++;
			p2++;
		}
	}
	return (*this);
}

template<typename Tp>
Matrix<Tp> Matrix<Tp>::operator + (const Matrix<Tp> &mat)
{
	try
	{
		iEqual(*this, mat);

		Matrix<Tp> add(mat.rows, mat.cols);
		Tp *p1 = pMat;
		Tp *p2 = mat.pMat;
		Tp *padd = add.pMat;

		for (int i = 0; i < mat.rows; i++)
		{
			for (int j = 0; j < mat.cols; j++)
			{
				*padd = *p1 + *p2;
				p1++;
				p2++;
				padd++;
			}
		}
		return add;
	}
	catch (noEqual)
	{
		cout << endl << "异常：不是同型矩阵，不能相加！" << endl;
		Matrix<Tp> add;
		return add;
	}
}

template<typename Tp>
Matrix<Tp> Matrix<Tp>::operator - (const Matrix<Tp> &mat)
{
	try
	{
		iEqual(*this, mat);

		Matrix<Tp> add(mat.rows, mat.cols);
		Tp *p1 = pMat;
		Tp *p2 = mat.pMat;
		Tp *padd = add.pMat;

		for (int i = 0; i < mat.rows; i++)
		{
			for (int j = 0; j < mat.cols; j++)
			{
				*padd = *p1 - *p2;
				p1++;
				p2++;
				padd++;
			}
		}
		return add;
	}
	catch (noEqual)
	{
		cout << endl << "异常：不是同型矩阵，不能相减！" << endl;
		Matrix<Tp> add;
		return add;
	}
}

template<typename Tp>
Matrix<Tp> Matrix<Tp>::operator * (const Matrix<Tp> &mat)
{
	try
	{
		iCanMultiply(*this, mat);

		Matrix<Tp> tar(rows, mat.cols);
		Tp *p1 = pMat;
		Tp *p2 = mat.pMat;
		Tp *ptar = tar.pMat;

		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < mat.cols; j++)
			{
				Tp sum_buf = 0;
				for (int k = 0; k < cols; k++)
				{
					sum_buf += p1[i*cols + k] * p2[k*mat.cols + j];
				}
				ptar[i*tar.cols + j] = sum_buf;
			}
		}
		return tar;
	}
	catch (cannotMultiply)
	{
		cout << endl << "异常：行数和列数的条件不满足，不能相乘！" << endl;
		Matrix<Tp> tar;
		return tar;
	}
}

template<typename Tp>
Matrix<Tp> Matrix<Tp>::operator / (const Tp num)
{
	Matrix<Tp> add(rows, cols);
	Tp *pmat = pMat;
	Tp *padd = add.pMat;

	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			*padd = *pmat / num;
			pmat++;
			padd++;
		}
	}
	return add;
}

template<typename Tp>
Matrix<Tp> Matrix<Tp>::operator * (const Tp num)
{
	Matrix<Tp> add(rows, cols);
	Tp *pmat = pMat;
	Tp *padd = add.pMat;

	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			*padd = *pmat * num;
			pmat++;
			padd++;
		}
	}
	return add;
}

template<typename Tp>
void Matrix<Tp>::operator () (const int row, int col)
{
	release();
	rows = row;
	cols = col;
	create();
}



template<typename Tp>
Matrix<Tp> operator - (const Matrix<Tp> &mat)
{
	Matrix<Tp> n_mat(mat.rows, mat.cols);
	Tp *pn = n_mat.pMat;
	Tp *p = mat.pMat;
	for (int i = 0; i < mat.rows; i++)
	{
		for (int j = 0; j < mat.cols; j++)
		{
			*pn = -*p;
			pn++;
			p++;
		}
	}
	return n_mat;
}

template<typename Tp>
Matrix<Tp> operator + (const Matrix<Tp> &mat)
{
	return mat;
}



template<typename Tp>
Matrix<Tp>::Matrix() : pMat(nullptr)
{
	rows = 0;
	cols = 0;
	//create();
}

template<typename Tp>
Matrix<Tp>::Matrix(const int &row, const int &col) : pMat(nullptr)
{
	rows = row;
	cols = col;
	create();
}

template<typename Tp>
void Matrix<Tp>::create()
{
	pMat = new Tp[rows*cols];
	Tp *p = pMat;
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			*p = 0.0;
			p++;
		}
	}
}

template<typename Tp>
void Matrix<Tp>::release()
{
	delete pMat;
}


template<typename Tp>
void Matrix<Tp>::print()
{
	cout << "[" << endl;
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			printf( "\t%.2f", pMat[i*cols + j]);
		}
		cout << endl;
	}
	cout << "]" << endl;
}

template<typename Tp>
Matrix<Tp> Matrix<Tp>::zeros()
{
	Tp *p = pMat;
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			*p = 0.0;
			p++;
		}
	}
	return (*this);
}

template<typename Tp>
Matrix<Tp> Matrix<Tp>::identity()
{
	Tp *p = pMat;
	if (cols != rows)
	{
		return (*this);
	}
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			if (i == j)
			{
				*p = 1.0;
			}
			else
			{
				*p = 0.0;
			}
			p++;
		}
	}
	return (*this);
}

template<typename Tp>
Matrix<Tp> Matrix<Tp>::ones()
{
	Tp *p = pMat;
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			*p = 1.0;
			p++;
		}
	}
	return (*this);
}

template<typename Tp>
Matrix<Tp> Matrix<Tp>::transpose()
{
	Matrix<Tp> tar(cols, rows);
	Tp *ptar = tar.pMat;
	Tp *p = pMat;
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			ptar[j*rows + i] = p[i*cols + j];
		}
	}
	return tar;
}


template<typename Tp>
void Matrix<Tp>::input()
{
	Tp *p = pMat;
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			cin >> (*p);
			p++;
		}
	}
}

template<typename Tp>
void Matrix<Tp>::iEqual(const Matrix<Tp> &mat1, const Matrix<Tp> &mat2)
{
	if (mat1.rows != mat2.rows || mat1.cols != mat2.cols)
	{
		throw noEqual();
	}
}

template<typename Tp>
void Matrix<Tp>::iCanMultiply(const Matrix<Tp> &mat1, const Matrix<Tp> &mat2)
{
	if (mat1.cols != mat2.rows)
	{
		throw cannotMultiply();
	}
}

/*
//	@名称：	RowOp
//	@参数：	op_type		-> 运算方式，有换法变换、倍法变换、消法变换，取值于枚举体RorCOpType_
//			row1		-> 第一行r1
//			row2		-> 第二行r2
//			k			-> 倍法变换下的系数
//	@说明：	换法变换—— r1与r2行互换	EXCHANGE_
//			倍法变换—— r1*k			MULTIPLY_
//			消法变换—— r1+r2*k		ELIMINATE_
*/
template<typename Tp>
Matrix<Tp> Matrix<Tp>::RowOp(RorCOpType_1 op_type, int row1, int row2)
{
	Tp *p1, *p2;
	p1 = p2 = pMat;
	if (row1 < 1)
	{
		row1 = 1;
	}
	else if (row1 > rows)
	{
		row1 = rows;
	}
	if (row2 < 1)
	{
		row2 = 1;
	}
	else if (row2 > rows)
	{
		row2 = rows;
	}
	row1 -= 1;
	row2 -= 1;
	p1 += row1 * cols;
	p2 += row2 * cols;

	Tp ele_buf;
	for (int i = 0; i < cols; i++)
	{
		ele_buf = *p1;
		*p1 = *p2;
		*p2 = ele_buf;
		p1++;
		p2++;
	}
	return (*this);
}
template<typename Tp>
Matrix<Tp>	Matrix<Tp>::RowOp(RorCOpType_2 op_type, int row1, double k)
{
	Tp *p1;
	p1 = pMat;
	if (row1 < 1)
	{
		row1 = 1;
	}
	else if (row1 > rows)
	{
		row1 = rows;
	}

	row1 -= 1;
	p1 += row1 * cols;

	for (int i = 0; i < cols; i++)
	{
		*p1 *= k;
		p1++;
	}
	return (*this);
}
template<typename Tp>
Matrix<Tp>	Matrix<Tp>::RowOp(RorCOpType_2 op_type, int row1, double k, int row2)
{
	Tp *p1, *p2;
	p1 = p2 = pMat;
	if (row1 < 1)
	{
		row1 = 1;
	}
	else if (row1 > rows)
	{
		row1 = rows;
	}
	if (row2 < 1)
	{
		row2 = 1;
	}
	else if (row2 > rows)
	{
		row2 = rows;
	}
	row1 -= 1;
	row2 -= 1;
	p1 += row1 * cols;
	p2 += row2 * cols;

	for (int i = 0; i < cols; i++)
	{
		*p1 += k * (*p2);
		p1++;
		p2++;
	}
	return (*this);
}

/*
//	@名称：	ColOp
//	@参数：	op_type		-> 运算方式，有换法变换、倍法变换、消法变换，取值于枚举体RorCOpType_
//			col1		-> 第一列c1
//			col2		-> 第二列c
//			k			-> 倍法变换下的系数
//	@说明：	换法变换—— c1与c2列互换	EXCHANGE_
//			倍法变换—— r2*k			MULTIPLY_
//			消法变换—— c1+c2*k		ELIMINATE_
*/
template<typename Tp>
Matrix<Tp> Matrix<Tp>::ColOp(RorCOpType_1 op_type, int col1, int col2)
{
	Tp *p1, *p2;
	p1 = p2 = pMat;
	if (col1 < 1)
	{
		col1 = 1;
	}
	else if (col1 > cols)
	{
		col1 = cols;
	}
	if (col2 < 1)
	{
		col2 = 1;
	}
	else if (col2 > cols)
	{
		col2 = cols;
	}
	col1 -= 1;
	col2 -= 1;
	p1 += col1;
	p2 += col2;

	Tp ele_buf;
	for (int i = 0; i < cols; i++)
	{
		ele_buf = *p1;
		*p1 = *p2;
		*p2 = ele_buf;
		p1 += rows;
		p2 += rows;
	}
	return (*this);
}
template<typename Tp>
Matrix<Tp>	Matrix<Tp>::ColOp(RorCOpType_2 op_type, int col1, double k)
{
	Tp *p1;
	p1 = pMat;
	if (col1 < 1)
	{
		col1 = 1;
	}
	else if (col1 > cols)
	{
		col1 = cols;
	}
	
	col1 -= 1;
	p1 += col1;

	for (int i = 0; i < cols; i++)
	{
		*p1 *= k;
		p1 += rows;
	}
	return (*this);
}
template<typename Tp>
Matrix<Tp>	Matrix<Tp>::ColOp(RorCOpType_2 op_type, int col1, double k, int col2)
{
	Tp *p1, *p2;
	p1 = p2 = pMat;
	if (col1 < 1)
	{
		col1 = 1;
	}
	else if (col1 > cols)
	{
		col1 = cols;
	}
	if (col2 < 1)
	{
		col2 = 1;
	}
	else if (col2 > cols)
	{
		col2 = cols;
	}
	col1 -= 1;
	col2 -= 1;
	p1 += col1;
	p2 += col2;

	for (int i = 0; i < cols; i++)
	{
		*p1 += k*(*p2);
		p1 += rows;
		p2 += rows;
	}
	return (*this);
}

template<typename Tp>
long double Matrix<Tp>::Determinant()
{
	Matrix<Tp> mat;
	vector<int> diagonal_pos(rows, 0);
	int deal_count = 0;
	long double det_val = 1;
	long double min;
	Tp *p = mat.pMat;
	mat = (*this);

	if (rows != cols)
	{
		return 0;
	}
	int i, j;
	for (j = 0; j < cols; j++)
	{
		//列元素不全为0，否则行列式为0
		for (i = 0, p = mat.pMat + j; i < rows; i++, p += rows)
		{
			if ((*p) != 0)
			{
				break;
			}
			else
			{
				return 0;
			}
		}
	}
	

	for ( j = 0; j < cols; j++)
	{
		int temp_row = deal_count;
		for (i = deal_count, p = mat.pMat + deal_count*cols +j; i < rows; i++, p += rows)
		{
			if (fabs(*p) > 1e-9)
			{
				min = *p;
			}
		}

		for (i = deal_count, p = mat.pMat + deal_count*cols + j; i < rows; i++, p += rows)
		{
			if (fabs(*p) > 1e-9 && fabs(*p) <= (fabs(min) + 1e-9))
			{
				min = *p;
				temp_row = i;
				diagonal_pos[i] = j;
			}
		}

		for ( i = deal_count, p = mat.pMat + deal_count*cols + j; i < rows; i++, p += rows)
		{
			if (i != temp_row && fabs(*p) > 1e-9)
			{
				long double k = (*p) / min;
				mat.RowOp(ELIMINATE_, i+1, -k, temp_row+1);
			}
		}
		mat.RowOp(EXCHANGE_, deal_count+1, temp_row+1);
		deal_count++;
	}
	
	p = mat.pMat;
	for (i = 0 ; i < rows; i++)
	{
		det_val *= p[i*cols + i];
	}

	return det_val;
}

template<typename Tp>
Matrix<Tp> Matrix<Tp>::Inverse()
{
	Matrix<Tp> mat;
	Matrix<Tp> inv(rows, cols);

	mat = (*this);
	inv.identity();

	vector<int> diagonal_pos(rows, 0);
	int deal_count = 0;
	long double det_val = 1;
	long double min;
	long double k;
	Tp *p = mat.pMat;

	if (rows != cols)
	{
		return (*this);
	}
	signed int i, j;

	//行列式为0，逆矩阵不存在
	if (abs(mat.Determinant()) < 1e-9)
	{
		return (*this);
	}

	for (j = 0; j < cols; j++)
	{
		int temp_row = deal_count;
		for (i = deal_count, p = mat.pMat + deal_count*cols + j; i < rows; i++, p += rows)
		{
			if (fabs(*p) > 1e-9)
			{
				min = *p;
			}
		}

		for (i = deal_count, p = mat.pMat + deal_count*cols + j; i < rows; i++, p += rows)
		{
			if (fabs(*p) > 1e-9 && fabs(*p) <= (fabs(min) + 1e-9))
			{
				min = *p;
				temp_row = i;
				diagonal_pos[i] = j;
			}
		}

		for (i = deal_count, p = mat.pMat + deal_count*cols + j; i < rows; i++, p += rows)
		{
			if (i != temp_row && fabs(*p) > 1e-9)
			{
				k = (*p) / min;
				mat.RowOp(ELIMINATE_, i + 1, -k, temp_row + 1);
				inv.RowOp(ELIMINATE_, i + 1, -k, temp_row + 1);
			}
		}
		mat.RowOp(EXCHANGE_, deal_count + 1, temp_row + 1);
		inv.RowOp(EXCHANGE_, deal_count + 1, temp_row + 1);
		deal_count++;
	}

	p = mat.pMat;
	for (j = cols - 1; j > 0; j--)
	{
		long double buf_val = p[j*cols + j];
		for (i = 0; i < j; i++)
		{
			k = p[cols * i + j] / buf_val;
			mat.RowOp(ELIMINATE_, i + 1, -k, j + 1);
			inv.RowOp(ELIMINATE_, i + 1, -k, j + 1);
		}
	}
	for (i = 0; i < rows; i++)
	{
		k = p[i*cols + i];
		inv.RowOp(MULTIPLY_, i + 1, 1.0 / k);
	}
	return inv;
	
}


#endif // ! _CHAPTER1_H_
