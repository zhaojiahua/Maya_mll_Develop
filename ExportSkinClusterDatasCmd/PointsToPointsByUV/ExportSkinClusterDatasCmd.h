#pragma once

#include <math.h>
#include <maya/MPxCommand.h>
#include <maya/MStatus.h>
#include <maya/MSelectionList.h>
#include <maya/MItSelectionList.h>
#include <maya/MFnMesh.h>
#include <maya/MFnDagNode.h>
#include <maya/MPlug.h>
#include <maya/MPlugArray.h>
#include <maya/MDagPath.h>
#include <maya/MDagPathArray.h>
#include <maya/MPoint.h>
#include <maya/MPointArray.h>
#include <maya/MVectorArray.h>
#include <maya/MIntArray.h>
#include <maya/MDoubleArray.h>
#include <maya/MObjectArray.h>
#include <maya/MArgList.h>
#include <maya/MArgParser.h>
#include <maya/MSyntax.h>
#include <maya/MItMeshPolygon.h>
#include <maya/MItDependencyNodes.h>
#include <maya/MItGeometry.h>
#include <maya/MFnSkinCluster.h>
#include <maya/MFnSingleIndexedComponent.h>
#include <maya/MIOStream.h>
#define CheckError(stat,msg)		\
			if(MS::kSuccess != stat) {	\
				displayError(msg);		\
				continue;						\
			}

class ExportSkinClusterDatas : public MPxCommand
{

public:
	ExportSkinClusterDatas();
	virtual		~ExportSkinClusterDatas();

	//指定dagpath,获取其skinCluster节点
	MStatus		getSkinCluster(const MObject& inNode, MObject& skinCluster);
	//解析参数
	MStatus		parseArgs(const MArgList& args);
	MStatus		doIt( const MArgList& args)override;
	MStatus		redoIt()override;
	MStatus		undoIt()override;
	bool		isUndoable() const override;

	static		void* creator();

private:
	// Store the data you will need to undo the command here
	MString filePathStr;//指定文件夹路径

};
