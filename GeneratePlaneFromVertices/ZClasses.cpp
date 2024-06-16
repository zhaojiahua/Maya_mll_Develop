#include "ZClasses.h"

Z5Matrix::Z5Matrix()
{
}

Z5Matrix::Z5Matrix(double* in25list)
{
}

Z5Matrix::Z5Matrix(Z5Vector* in5vects)
{
}

void Z5Matrix::SetElementsAll(double* in25list)
{
}

void Z5Matrix::SetElement(int index_r, int index_c, double value)
{
}

void Z5Matrix::SetElementsByZ5Vector(int index, Z5Vector in5vec)
{
}

Z5Vector Z5Matrix::GetElementsByZ5Vector(int index)
{
	return Z5Vector();
}

double Z5Matrix::GetElement(int index_r, int index_c)
{
	return 0.0;
}

void Z5Matrix::MakeIdentity()
{
}

Z5Matrix Z5Matrix::Transpose()
{
	return Z5Matrix();
}

double* Z5Matrix::SubMatrix(double* indata, int index_r, int index_c)
{
	return nullptr;
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
	return Z5Matrix();
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
}

Z5Vector::Z5Vector()
{
	for (int i = 0; i < 5; ++i) {	mateData[i] = 0;}
}

Z5Vector::Z5Vector(double* in5list)
{
	for (int i = 0; i < 5; ++i) { mateData[i] = in5list[i]; }
}

void Z5Vector::SetElements(double* in5list)
{
	for (int i = 0; i < 5; ++i) { mateData[i] = in5list[i]; }
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
