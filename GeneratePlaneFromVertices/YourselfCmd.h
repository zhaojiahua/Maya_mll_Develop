#pragma once
#include <maya/MPxCommand.h>
#include <maya/MGlobal.h>
#include <maya/MSelectionList.h>
#include <maya/MItSelectionList.h>
#include <maya/MStringArray.h>
#include <maya/MString.h>
#include <maya/MDagPath.h>
#include <maya/MDagPathArray.h>
#include <maya/MArgList.h>
#include <maya/MObject.h>
#include <maya/MFnDagNode.h>
#include <maya/MFnMesh.h>
#include <maya/MPoint.h>
#include <maya/MPointArray.h>
#include <maya/MSyntax.h>
#include <maya/MArgParser.h>
#include <maya/MFnNurbsCurve.h>
#include <maya/MFnTransform.h>
#include <maya/MColor.h>
#include <maya/MColorArray.h>
#include <maya/MIntArray.h>
#include <maya/MVector.h>
#include <maya/MVectorArray.h>
#include <maya/MQuaternion.h>
#include <maya/MMatrix.h>
#include <maya/MPlugArray.h>
#include <maya/MPlug.h>
#include <maya/MMaterial.h>
#include <maya/MHardwareRenderer.h>
#include <maya/MFragmentManager.h>
#include <maya/MFnNurbsCurveData.h>
#include <maya/MFnDependencyNode.h>

#include "ZClasses.h"

//Maya��������
class YourselfCommand :public MPxCommand
{
public:
	YourselfCommand();
	virtual		~YourselfCommand();

	MStatus		doIt(const MArgList&);
	MStatus		redoIt();
	MStatus		undoIt();
	bool		isUndoable() const;

	static		void* creator();

private:// Store the data you will need to undo the command here
	
	MSelectionList selections;
	unsigned int crvptnum = 10;
	unsigned int rowcurveSpans = 5;
	double bextend = 0.05;

	MDagPathArray resultsDagpaths;

	//��mesh�ŵ���������,�����ŵ����ʵĴ�С(�������������е�mesh��,mesh������λ��,�����ű���,˳�����չ�����������)
	MStatus FitMesh(const MPointArray& inpoints,bool makeXsp, MPointArray&, MPoint&, double&, unsigned int&, MDoubleArray&);
	//������ѡ���Mesh,���᷵�ظ�Mesh������������(0��x��,1��y��,2��z��),��ӳ�䵽��Ӧ������ƽ���϶�ά��,�Լ��ڸ������ϵĿ��(���ֵ����Сֵ֮��)
	MStatus Get2DPoints(const MPointArray& allpoints, const unsigned int& saxies, MDoubleArray& oarray1, MDoubleArray& oarray2);
	//����һ��ֵ,�����ֵ���մ�С���������˳�򰲲嵽��Ӧ��λ��
	MStatus SortedAppend(double invalue, MDoubleArray& darray);
	//����һ�Ѷ�ά���ȡ���������
	MStatus GetAugmentedMatrix(MDoubleArray inpts, Z5Matrix& oMat, Z5Vector& oVec);
	//��һ��MDoubleArray�۰��ֳ�����MDoubleArray
	MStatus HalfSplit(MDoubleArray inpts12, MDoubleArray& inpts1, MDoubleArray& inpts2);
	//�������ʽϵ��(����ƽ��ռ�),��չ����,��չ����ȡֵ��Χ,�ͷֱ���,�������ߵ�EP������
	MPointArray GenerateCrvEPs(Z5Vector inc1, Z5Vector inc2, unsigned int spreadAxies, const MDoubleArray& mmv, unsigned int epnum,double inbe);
	
	//����㼯��һ������,��Щ�㼯�ᱻ�����������ߵķ���ָ�����������Ӧ�Ķ�����Ӧ(���ص����������)
	MIntArray* SlitedPointsByCrv(const MPointArray& inpoints, const MFnNurbsCurve& incrv);
	//У��Normals(��˳һ�������ϵ�����Normals,ʹ���ǳ���һ��)(Ҫ������֮һһ��Ӧ������)
	MStatus StraightenNormals(MVectorArray& innormals, const MVectorArray& intangents);
	//����㼯�Ϻ�Ҫ���ɵ�EP�����������Ӧ������������ߵ�EP��
	MPointArray GenerateMatchedCurveEPs(const MPointArray& allpoints, bool makeXsp = false, unsigned int epcount = 5);
	//����EP������������������,�������ɵ�����
	MDagPath GenerateEPCurve(const MPointArray& ineps, MString inname, MPointArray& unieps);
	//�������EP��,�������polygon�Ĳ���
	ZGenMeshParam GenMeshParam(MPointArray* inpointArrays, unsigned int ptarrayNum);
	//���봴��polygon�Ĳ���������Ӧ��polygon
	MDagPath ZGenMesh(const ZGenMeshParam& inmeshparam, MString inname);
	//����EP����,����Щ������ƽ��(˳��ѵ���0���������ȥ)
	MPointArray* SmoothUniEPS(MPointArray* inunieps, unsigned int& inptarrayNum);
};