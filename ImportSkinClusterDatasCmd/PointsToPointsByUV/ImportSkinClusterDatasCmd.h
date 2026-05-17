#pragma once

#include <math.h>
#include <maya/MPxCommand.h>
#include <maya/MStatus.h>
#include <maya/MArgDatabase.h>
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
#include <maya/MItDependencyGraph.h>
#include <maya/MItGeometry.h>
#include <maya/MFnSkinCluster.h>
#include <maya/MFnSingleIndexedComponent.h>
#include <maya/MIOStream.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#define CheckError(stat,msg)		\
			if(MS::kSuccess != stat) {	\
				displayError(msg);		\
				continue;						\
			}

using namespace std;
class ImportSkinClusterDatas : public MPxCommand
{
public:
	ImportSkinClusterDatas();
	virtual		~ImportSkinClusterDatas();

	MStatus		doIt( const MArgList& args)override;
	MStatus		redoIt()override;
	MStatus		undoIt()override;
	bool		isUndoable() const override;

	static		void* creator();

	static		MSyntax newSyntax();

private:
	// Store the data you will need to undo the command here
	MString m_filePathStr;//指定csv文件夹路径
	MSelectionList m_selectionlist;
	MObject m_skinClusterNode{ MObject::kNullObj };
	MObject m_skinClusterNode_old{ MObject::kNullObj };
	MDagPath m_meshDagPath;
	MObject m_meshComponent;
	MStringArray m_jointsNames;
	vector<MDoubleArray> m_vertexWeights;		//存储每个顶点的新的权重数据,用于重做
	vector<MDoubleArray> m_vertexWeights_old;		//存储每个顶点旧的权重数据,用于撤销操作

	//指定DagNode,获取其skinCluster节点
	MStatus		getSkinCluster(MObject inNode, MObject& skinCluster);
	//解析参数
	MStatus		parseArgs(const MArgList& args);
	//读取.CSV文件,将其首行的joints名称存储在m_jointsNames,其对应的的顶点权重数据存储在m_vertexWeights
	MStatus		readCSVFile(const MString& fullpath);
	//创建蒙皮
	MStatus		createSkinCluster();
	//应用权重数据
	MStatus		applyWeights(MObject& inSkinCluster, const vector<MDoubleArray>& inVertexWeights, vector<MDoubleArray>& inVertexWeights_old);
};
