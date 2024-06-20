#include "ZClasses.h"

Z5Matrix::Z5Matrix()
{
	mateData = MDoubleArray(25, 0.0);
	for (int i = 0; i < 5; ++i) {
		for (int j = 0; j < 5; ++j) {
			if (i == j)mateData[i * 5 + j] = 1;
		}
	}
}

Z5Matrix::Z5Matrix(double* in25list)
{
	mateData = MDoubleArray(in25list, 25);
}

Z5Matrix::Z5Matrix(MDoubleArray in25array)
{
	mateData = in25array;
}

Z5Matrix::Z5Matrix(Z5Vector* in5vects)
{
	mateData = MDoubleArray(25, 0.0);
	for (int i = 0; i < 5; ++i) {
		for (int j = 0; j < 5; ++j) {
			mateData[i * 5 + j] = in5vects[i].GetElement(j);
		}
	}
}

Z5Matrix::~Z5Matrix()
{}

void Z5Matrix::SetElementsAll(double* in25list)
{
	mateData = MDoubleArray(in25list, 25);
}
void Z5Matrix::SetElementsAll(MDoubleArray in25array)
{
	mateData = in25array;
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

MDoubleArray Z5Matrix::GetMateData()
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

MDoubleArray Z5Matrix::SubMatrix(MDoubleArray inmat, int index_r, int index_c)
{
	int matsize = sqrt(inmat.length());
	MDoubleArray newmate;
	if (matsize < 2) { MGlobal::displayError("matrix size < 2"); return MDoubleArray(); }
	else {
		for (int i = 0; i < matsize; ++i) {
			for (int j = 0; j < matsize; ++j) {
				if (i != index_r && j != index_c)newmate.append(inmat[i * matsize + j]);
			}
		}
		return newmate;
	}
}

Z5Matrix Z5Matrix::operator*(Z5Matrix in5mat)
{
	Z5Matrix newmat;
	for (int i = 0; i < 5; ++i) {
		for (int j = 0; j < 5; ++j) {
			double tresult = 0;
			for (int k = 0; k < 5; ++k) {
				tresult += (GetElement(i, k) * in5mat.GetElement(k, j));
			}
			newmat.SetElement(i, j, tresult);
		}
	}
	return newmat;
}

Z5Vector Z5Matrix::operator*(Z5Vector in5vec)
{
	Z5Vector tvector;
	for (int r = 0; r < 5; ++r) {
		double tvalue = 0;
		for (int i = 0; i < 5; ++i) {
			tvalue += (GetElement(r, i) * in5vec.GetElement(i));
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

double Z5Matrix::Determinant()
{
	return Z5Matrix::Determinant(mateData);
}

double Z5Matrix::Determinant(MDoubleArray inmate)
{
	int matsize = sqrt(inmate.length());
	if (matsize == 2)return inmate[0] * inmate[3] - inmate[1] * inmate[2];
	else {
		double tvalue = 0;
		for (int i = 0; i < matsize; ++i) {
			tvalue += (inmate[i] * Z5Matrix::Determinant(Z5Matrix::SubMatrix(inmate, 0, i)) * pow(-1, i));
		}
		return tvalue;
	}
}

Z5Matrix Z5Matrix::AdjointMatrix()
{
	Z5Matrix tmat;
	for (int i = 0; i < 5; ++i) {
		for (int j = 0; j < 5; ++j) {
			double dsyzs = pow(-1, i + j) * Z5Matrix::Determinant(Z5Matrix::SubMatrix(mateData, i, j));
			tmat.SetElement(j, i, dsyzs);
		}
	}
	return tmat;
}

Z5Matrix Z5Matrix::InverseMatrix()
{
	return 1.0 / Determinant() * AdjointMatrix();
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
	MGlobal::displayInfo("");
}

void Z5Matrix::PrintMateDatas(MDoubleArray indatas)
{
	int insize = sqrt(indatas.length());
	for (int i = 0; i < insize; ++i) {
		MString tstr = MString("[");
		for (int j = 0; j < insize; ++j) {
			if (j != 0)tstr += ",";
			tstr += indatas[i * insize + j];
		}
		tstr += "]";
		MGlobal::displayInfo(tstr);
	}
	MGlobal::displayInfo("");
}

void Z5Matrix::PrintMateDatas(double indata)
{
	MGlobal::displayInfo(MString("") + indata);
	MGlobal::displayInfo("");
}

Z5Vector::Z5Vector()
{
	mateData = MDoubleArray(5, 0.0);
}

Z5Vector::Z5Vector(double* in5list)
{
	mateData = MDoubleArray(in5list, 5);
}

Z5Vector::Z5Vector(MDoubleArray inarry)
{
	mateData = inarry;
}

void Z5Vector::SetElements(double* in5list)
{
	mateData = MDoubleArray(in5list, 5);
}

void Z5Vector::SetElements(MDoubleArray inarry)
{
	mateData = inarry;
}

void Z5Vector::SetElement(int index, double invalue)
{
	mateData[index] = invalue;
}

double Z5Vector::GetElement(int index)
{
	return mateData[index];
}

MDoubleArray Z5Vector::GetMateData()
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

Z5Matrix operator*(double ind, Z5Matrix inm)
{
	return inm * ind;
}

MDoubleArray operator+(MDoubleArray a1, MDoubleArray a2)
{
	MDoubleArray  nea = MDoubleArray();
	for (int i = 0; i < a1.length(); ++i)nea.append(a1[i]);
	for (int i = 0; i < a2.length(); ++i)nea.append(a2[i]);
	return nea;
}
