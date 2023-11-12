#pragma once
#include <maya/MPxCommand.h>
#include <maya/MSelectionList.h>
#include <vector>
#include <maya/MDagPathArray.h>

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

private:
	//for debug
	MString tempStr;
	MSelectionList selectionList;
	//根据给定的根节点,遍历获取其下的所有子节点(如果某个根节点下有多个子节点,那么只有放在第一位的子节点才生效)
	MStatus GetAllChildrenPosition(const MDagPath& rootPath, MVectorArray& outArray);
	std::vector<MVectorArray> rootPoss;
	MDagPathArray jntsPaths;
	MObject outMesh;
};