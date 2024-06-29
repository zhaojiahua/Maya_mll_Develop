#pragma once
#include <maya/MGlobal.h>
#include <maya/MDoubleArray.h>
#include <maya/MRandom.h>
#include <maya/MPointArray.h>
#include <maya/MPoint.h>
#include <maya/MIntArray.h>

using namespace std;

class Z5Matrix;

class Z5Vector
{
	MDoubleArray mateData;
public:
	Z5Vector();
	Z5Vector(double* in5list);
	Z5Vector(MDoubleArray inarry);

	void SetElements(double* in5list);
	void SetElements(MDoubleArray inarry);
	void SetElement(int index, double invalue);
	double GetElement(int index);
	MDoubleArray GetMateData();

	Z5Vector operator+ (Z5Vector inVec);
	Z5Vector operator* (double inS);
	Z5Vector operator* (Z5Vector inVec);//非向量乘法
	Z5Vector operator* (Z5Matrix inMat);

	//打印函数(方便Debug)
	void Print();
};

class Z5Matrix
{
	MDoubleArray mateData;
public:
	Z5Matrix();
	Z5Matrix(double* in25list);
	Z5Matrix(MDoubleArray in25array);
	Z5Matrix(Z5Vector* in5vects);

	~Z5Matrix();

	void SetElementsAll(double* in25list);
	void SetElementsAll(MDoubleArray in25array);
	void SetElement(int index_r, int index_c, double value);
	void SetElementsByZ5Vector(int index, Z5Vector in5vec);
	Z5Vector GetElementsByZ5Vector(int index);
	MDoubleArray GetMateData();
	double GetElement(int index_r, int index_c);
	//使该矩阵变成单位矩阵
	void MakeIdentity();
	//返回该矩阵的转置矩阵
	Z5Matrix Transpose();
	//返回指定元素坐标的子矩阵(去掉该行该列后得到的matsize-1阶矩阵)(返回的是mateData数组)(要传入原始数据和方阵的大小)
	static MDoubleArray SubMatrix(MDoubleArray inmat, int index_r, int index_c);

	Z5Matrix operator*(Z5Matrix in5mat);
	Z5Vector operator*(Z5Vector in5vec);
	Z5Matrix operator*(double inval);

	//返回此矩阵的行列式的值
	double Determinant();
	//输入方阵返回该矩阵的行列式的值
	static double Determinant(MDoubleArray inmate);
	//返回该矩阵的伴随矩阵
	Z5Matrix AdjointMatrix();
	//返回该矩阵的逆矩阵
	Z5Matrix InverseMatrix();

	//打印函数(方便Debug)
	void Print();
	//打印matedata数组(需要指定方阵的大小)
	static void PrintMateDatas(MDoubleArray indatas);
	static void PrintMateDatas(double indata);
};

Z5Matrix operator*(double, Z5Matrix);

MDoubleArray operator+(MDoubleArray a1, MDoubleArray a2);

struct ZGenMeshParam
{
	int numVertices;
	int numPolygons;
	MPointArray vertexArray;
	MIntArray polygonCounts;
	MIntArray polygonConnects;
	MObject parentOrOwner;
	MStatus* ReturnStatus;

	ZGenMeshParam();
};