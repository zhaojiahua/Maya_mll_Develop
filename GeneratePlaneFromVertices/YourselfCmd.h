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

	//传入所选择的Mesh,它会返回该Mesh所在延伸轴向(0是x轴,1是y轴,2是z轴),和映射到对应的两个平面上二维点,以及在该轴向上的跨度(最大值与最小值之差)
	MStatus Get2DPoints(const MDagPath& inobj, int&, MDoubleArray&, double&);

};