#pragma once
#include <maya/MPxCommand.h>
#include <maya/MGlobal.h>
#include <maya/MSelectionList.h>
#include <maya/MItSelectionList.h>
#include <maya/MStringArray.h>
#include <maya/MString.h>
#include <maya/MDagPath.h>
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
	double bextend = 0.05;

	//传入所选择的Mesh,它会返回该Mesh所在延伸轴向(0是x轴,1是y轴,2是z轴),和映射到对应的两个平面上二维点,以及在该轴向上的跨度(最大值与最小值之差)
	MStatus Get2DPoints(const MDagPath& inobj, unsigned int&, MDoubleArray&, MDoubleArray&);
	//传入一个值,把这个值按照从小到大的排列顺序安插到相应的位置
	MStatus SortedAppend(double invalue, MDoubleArray& darray);
	//输入一堆二维点获取其曾广矩阵
	MStatus GetAugmentedMatrix(MDoubleArray inpts, Z5Matrix& oMat, Z5Vector& oVec);
	//把一个MDoubleArray折半拆分成两个MDoubleArray
	MStatus HalfSplit(MDoubleArray inpts12, MDoubleArray& inpts1, MDoubleArray& inpts2);
	//传入多项式系数(两个平面空间),延展轴向,延展轴向取值范围,和分辨率,返回曲线的EP点坐标
	MStatus GenerateCrvEPs(Z5Vector inc1, Z5Vector inc2, unsigned int spreadAxies, const MDoubleArray& mmv, unsigned int epnum, MPointArray& oPoints);
};