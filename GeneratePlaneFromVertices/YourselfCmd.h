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

//Maya的命令类
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

	//把mesh放到世界中心,且缩放到合适的大小(输出调整后的所有的mesh点,mesh的中心位置,和缩放比例,顺便把延展轴向给返回了)
	MStatus FitMesh(const MPointArray& inpoints,bool makeXsp, MPointArray&, MPoint&, double&, unsigned int&, MDoubleArray&);
	//传入所选择的Mesh,它会返回该Mesh所在延伸轴向(0是x轴,1是y轴,2是z轴),和映射到对应的两个平面上二维点,以及在该轴向上的跨度(最大值与最小值之差)
	MStatus Get2DPoints(const MPointArray& allpoints, const unsigned int& saxies, MDoubleArray& oarray1, MDoubleArray& oarray2);
	//传入一个值,把这个值按照从小到大的排列顺序安插到相应的位置
	MStatus SortedAppend(double invalue, MDoubleArray& darray);
	//输入一堆二维点获取其曾广矩阵
	MStatus GetAugmentedMatrix(MDoubleArray inpts, Z5Matrix& oMat, Z5Vector& oVec);
	//把一个MDoubleArray折半拆分成两个MDoubleArray
	MStatus HalfSplit(MDoubleArray inpts12, MDoubleArray& inpts1, MDoubleArray& inpts2);
	//传入多项式系数(两个平面空间),延展轴向,延展轴向取值范围,和分辨率,返回曲线的EP点坐标
	MPointArray GenerateCrvEPs(Z5Vector inc1, Z5Vector inc2, unsigned int spreadAxies, const MDoubleArray& mmv, unsigned int epnum,double inbe);
	
	//传入点集和一条曲线,这些点集会被沿着这条曲线的方向分割成与曲线相对应的段数对应(返回点的索引集合)
	MIntArray* SlitedPointsByCrv(const MPointArray& inpoints, const MFnNurbsCurve& incrv);
	//校正Normals(理顺一个曲线上的所有Normals,使他们朝向一致)(要传入与之一一对应的切线)
	MStatus StraightenNormals(MVectorArray& innormals, const MVectorArray& intangents);
	//传入点集合和要生成的EP点个数生成相应个数的拟合曲线的EP点
	MPointArray GenerateMatchedCurveEPs(const MPointArray& allpoints, bool makeXsp = false, unsigned int epcount = 5);
	//传入EP点和名称生成拟合曲线,返回生成的曲线
	MDagPath GenerateEPCurve(const MPointArray& ineps, MString inname, MPointArray& unieps);
	//传入多组EP点,输出创造polygon的参数
	ZGenMeshParam GenMeshParam(MPointArray* inpointArrays, unsigned int ptarrayNum);
	//传入创造polygon的参数生成相应的polygon
	MDagPath ZGenMesh(const ZGenMeshParam& inmeshparam, MString inname);
	//传入EP点阵,将这些点整理并平滑(顺便把等于0的曲线组除去)
	MPointArray* SmoothUniEPS(MPointArray* inunieps, unsigned int& inptarrayNum);
};