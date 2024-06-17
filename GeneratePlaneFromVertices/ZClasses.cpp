#include "ZClasses.h"
#include <math.h>

Z5Matrix::Z5Matrix()
{
	mateData = new double[25];
	for (int i = 0; i < 5; ++i) {
		for (int j = 0; j < 5; ++j) {
			if (i == j)mateData[i * 5 + j] = 1;
			else mateData[i * 5 + j] = 0;
		}
	}
}

Z5Matrix::Z5Matrix(double* in25list)
{
	mateData = in25list;
}

Z5Matrix::Z5Matrix(Z5Vector* in5vects)
{
	mateData = new double[25];
	for (int i = 0; i < 5; ++i) {
		for (int j = 0; j < 5; ++j) {
			mateData[i * 5 + j] = in5vects[i].GetElement(j);
		}
	}
}

Z5Matrix::~Z5Matrix()
{
	delete[] mateData;
}

void Z5Matrix::SetElementsAll(double* in25list)
{
	mateData = in25list;
}

void Z5Matrix::SetElement(int index_r, int index_c, double value)
{
	mateData[index_r * 5 + index_c] = value;
}

void Z5Matrix::SetElementsByZ5Vector(int index, Z5Vector in5vec)
{
	for (int i = 0; i < 5; ++i) {
		mateData[index * 5 + i] = in5vec.GetElement(i);
	}
}

Z5Vector Z5Matrix::GetElementsByZ5Vector(int index)
{
	double tempvalues[5];
	for (int i = 0; i < 5; ++i) {
		tempvalues[i] = mateData[index * 5 + i];
	}
	return Z5Vector(tempvalues);
}

double* Z5Matrix::GetMateData()
{
	return mateData;
}

double Z5Matrix::GetElement(int index_r, int index_c)
{
	return mateData[index_r * 5 + index_c];
}

void Z5Matrix::MakeIdentity()
{
	for (int i = 0; i < 5; ++i) {
		for (int j = 0; j < 5; ++j) {
			if (i == j)mateData[i * 5 + j] = 1;
			else mateData[i * 5 + j] = 0;
		}
	}
}

Z5Matrix Z5Matrix::Transpose()
{
	double tmate[25];
	for (int i = 0; i < 5; ++i) {
		for (int j = 0; j < 5; ++j) {
			tmate[i * 5 + j] = mateData[j * 5 + i];
		}
	}
	return Z5Matrix(tmate);
}

double* Z5Matrix::SubMatrix(double* indata, int matsize, int index_r, int index_c)
{
	if (matsize < 2) { MGlobal::displayError("matrix size < 2"); return NULL; }
	else {
		int tsize = (matsize - 1) * (matsize - 1);
		double* tmat = new double(tsize);
		for (int i = 0; i < index_r; ++i) {
			for (int j = 0; j < index_c; ++j) {
				tmat[i * (matsize - 1) + j] = indata[i * matsize + j];
			}
		}
		for (int i = index_r+1; i < matsize; ++i) {
			for (int j = 0; j < index_c; ++j) {
				tmat[(i-1) * (matsize - 1) + j] = indata[i * matsize + j];
			}
		}
		/*for (int i = 0; i < index_r; ++i) {
			for (int j = index_c+1; j < matsize; ++j) {
				tmat[i * (matsize - 1) + j - 1] = indata[i * matsize + j];
			}
		}
		for (int i = index_r + 1; i < matsize; ++i) {
			for (int j = index_c + 1; j < matsize; ++j) {
				tmat[(i - 1) * (matsize - 1) + j - 1] = indata[i * matsize + j];
			}
		}*/
		return tmat;
	}
}

Z5Matrix Z5Matrix::operator*(Z5Matrix in5mat)
{
	return Z5Matrix();
}

Z5Vector Z5Matrix::operator*(Z5Vector in5vec)
{
	Z5Vector tvector;
	for (int r = 0; r < 5; ++r) {
		double tvalue = 0;
		for (int i = 0; i < 5; ++i) {
			tvalue += (mateData[r * 5 + i] * in5vec.GetMateData()[i]);
		}
		tvector.SetElement(r, tvalue);
	}
	return tvector;
}

Z5Matrix Z5Matrix::operator*(double inval)
{
	Z5Matrix tmat = Z5Matrix();
	for (int i = 0; i < 5; ++i) {
		for (int j = 0; j < 5; ++j) {
			tmat.SetElement(i, j, GetElement(i, j) * inval);
		}
	}
	return tmat;
}

double Z5Matrix::M2X2_Det(double* invalues)
{
	return 0.0;
}

double Z5Matrix::M3X3_Det(double* invalues)
{
	return 0.0;
}

double Z5Matrix::M4X4_Det(double* invalues)
{
	return 0.0;
}

double Z5Matrix::Determinant()
{
	return 0.0;
}

Z5Matrix Z5Matrix::AdjointMatrix()
{
	return Z5Matrix();
}

Z5Matrix Z5Matrix::InverseMatrix()
{
	return Z5Matrix();
}

void Z5Matrix::Print()
{
	for (int i = 0; i < 5; ++i) {
		MString tstr = MString("[");
		for (int j = 0; j < 5; ++j) {
			if (j != 0)tstr += ",";
			tstr += mateData[i * 5 + j];
		}
		tstr += "]";
		MGlobal::displayInfo(tstr);
	}
}

void Z5Matrix::PrintMateDatas(double* indatas, int insize)
{
	for (int i = 0; i < insize; ++i) {
		MString tstr = MString("[");
		for (int j = 0; j < insize; ++j) {
			if (j != 0)tstr += ",";
			tstr += indatas[i * insize + j];
		}
		tstr += "]";
		MGlobal::displayInfo(tstr);
	}
}

Z5Vector::Z5Vector()
{
	for (int i = 0; i < 5; ++i) {	mateData[i] = 0;}
}

Z5Vector::Z5Vector(double* in5list)
{
	mateData = in5list;
}

void Z5Vector::SetElements(double* in5list)
{
	mateData = in5list;
}

void Z5Vector::SetElement(int index, double invalue)
{
	mateData[index] = invalue;
}

double Z5Vector::GetElement(int index)
{
	return mateData[index];
}

double* Z5Vector::GetMateData()
{
	return mateData;
}

Z5Vector Z5Vector::operator+(Z5Vector inVec)
{
	Z5Vector newVector;
	for (int i = 0; i < 5; ++i) {
		newVector.mateData[i] = mateData[i] + inVec.mateData[i];
	}
	return newVector;
}

Z5Vector Z5Vector::operator*(double inS)
{
	Z5Vector newVector;
	for (int i = 0; i < 5; ++i) {
		newVector.mateData[i] = mateData[i] * inS;
	}
	return newVector;
}

Z5Vector Z5Vector::operator*(Z5Vector inVec)
{
	Z5Vector newVector;
	for (int i = 0; i < 5; ++i) {
		newVector.mateData[i] = mateData[i] * inVec.mateData[i];
	}
	return newVector;
}

Z5Vector Z5Vector::operator*(Z5Matrix inMat)
{
	return inMat * Z5Vector(mateData);
}

void Z5Vector::Print()
{
	MString tstr("[");
	MGlobal::displayInfo(tstr + mateData[0] + "," + mateData[1] + "," + mateData[2] + "," + mateData[3] + "," + mateData[4] + "]");
}

