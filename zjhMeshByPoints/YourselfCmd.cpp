#include "YourselfCmd.h"
#include <maya/MGlobal.h>
#include <maya/MFnMesh.h>
#include <maya/MFloatPoint.h>
#include <maya/MFloatPointArray.h>
#include <maya/MIntArray.h>
#include <maya/MItSelectionList.h>
#include <maya/MDagPath.h>
#include <maya/MFnDagNode.h>
#include <maya/MFnTransform.h>
#include <maya/MFnSkinCluster.h>

YourselfCommand::YourselfCommand()
{}
YourselfCommand::~YourselfCommand()
{}

MStatus YourselfCommand::doIt(const MArgList&)
{
	MStatus stat = MS::kSuccess;
	MGlobal::getActiveSelectionList(selectionList);
	for (unsigned int i = 0; i < selectionList.length(); i++)
	{
		MDagPath tempath;
		selectionList.getDagPath(i, tempath);		//获取DagPath
		jntsPaths.append(tempath);
		MVectorArray tempVetorArray;
		MStringArray tempJointsname;
		stat = GetAllChildrenPosition(tempath, tempVetorArray, tempJointsname);
		rootPoss.push_back(tempVetorArray);

		for (int i=0;i< tempJointsname.length();++i)
		{
			jointsname.append(tempJointsname[i]);
		}
	}
	MGlobal::displayInfo(tempStr + "rootPoss.size: " + rootPoss.size());
	if (rootPoss.size() > 1) {
		for (int i = 1; i < rootPoss.size(); ++i) {
			if (rootPoss[i].length() != rootPoss[i - 1].length()) {
				MGlobal::displayWarning("each and every root joints children must be equal");//每个根骨骼的子骨骼必须数目相等
				break;
			}
		}
	}
	else {
		MGlobal::displayWarning("the root number at least 2 joints");
	}
	return redoIt();
}

MStatus YourselfCommand::redoIt()
{
	MStatus stat = MS::kSuccess;

	MFnMesh tempFnMesh;
	//int ptNum, faceNum;
	MFloatPointArray ptArrya;
	MIntArray faceCounts, faceConnects;
	//stat = GetParametersForCreator(rootPoss);
	unsigned int rootnumb = rootPoss.size();
	if (rootnumb > 1) {
		for (unsigned int i = 0; i < rootnumb; i++) {
			unsigned int next_i = (i + 1) % rootnumb;
			for (unsigned int j = 0; j < rootPoss[i].length(); j++) {
				ptArrya.append(MFloatPoint(rootPoss[i][j]));
				if (j < rootPoss[i].length() - 1 && j < rootPoss[next_i].length() - 1 && i < rootnumb - 1) {
					faceCounts.append(4);//每个面有4个点
					//4个点的序号
					faceConnects.append(i * rootPoss[i].length() + j);//polygon pt1
					faceConnects.append(i * rootPoss[i].length() + j + 1);//polygon pt2
					faceConnects.append(next_i * rootPoss[next_i].length() + j + 1);//polygon pt3
					faceConnects.append(next_i * rootPoss[next_i].length() + j);//polygon pt4
				}
			}
		}
		outMesh = tempFnMesh.create(ptArrya.length(), faceCounts.length(), ptArrya, faceCounts, faceConnects, MObject::kNullObj, &stat);
		MFnDagNode outMeshNode(outMesh);
		outMesh_name = outMeshNode.name();
		jointsname.append(outMesh_name);
	}
	else {
		MGlobal::displayWarning("the root number at least 2 joints");
	}
	MGlobal::selectByName(outMesh_name, MGlobal::kReplaceList);
	setResult(jointsname);
	return stat;
}

MStatus YourselfCommand::GetAllChildrenPosition(const MDagPath& rootPath, MVectorArray& outArray, MStringArray& outJointsnames)
{

	MFnTransform fnTransform(rootPath);	//用来获取当前节点的世界坐标
	outArray.append(fnTransform.getTranslation(MSpace::kWorld));
	MFnDagNode fnDagNode(rootPath);		//用来获取当前节点的子节点
	outJointsnames.append(fnDagNode.name());
	if (fnDagNode.childCount() == 0) return MS::kSuccess;
	else {
		MObject temp_child = fnDagNode.child(0);
		MFnDagNode child_fnDagNode(temp_child);
		MDagPath childPath;
		child_fnDagNode.getPath(childPath);
		GetAllChildrenPosition(childPath, outArray, outJointsnames);
	}
	return MS::kSuccess;
}

MStatus YourselfCommand::undoIt()
{
	MGlobal::deleteNode(outMesh);
	MGlobal::displayInfo("zjhMeshByPoints command undone!\n");
	return MS::kSuccess;
}

bool YourselfCommand::isUndoable() const
{

	return true;//可撤销的命令
}

void* YourselfCommand::creator()
{
	return new YourselfCommand();
}
