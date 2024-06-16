#include "YourselfCmd.h"


YourselfCommand::YourselfCommand()
{}
YourselfCommand::~YourselfCommand()
{}

MStatus YourselfCommand::doIt(const MArgList& arglist)
{
	MStatus stat = MS::kSuccess;
	MSelectionList allselections;
	MGlobal::getActiveSelectionList(allselections);
	MStatus suc;
	for (int i=0; i< allselections.length();++i) {
		MDagPath tdag;
		allselections.getDagPath(i, tdag);
		MFnDagNode tfnnode(tdag);
		if (tfnnode.child(0).apiType() == MFn::kMesh) {
			selections.add(tdag);
		}
	}

	return redoIt();
}

MStatus YourselfCommand::redoIt()
{
	MStatus suc;
	MItSelectionList itselections(selections);
	for (; itselections.isDone() != true; itselections.next()) {
		MDagPath tdag;
		itselections.getDagPath(tdag);
		MFnMesh tfnmesh(tdag);
		MPointArray tpoints;
		suc=tfnmesh.getPoints(tpoints, MSpace::kWorld);
		for (int i = 0; i < tpoints.length(); ++i){
			MGlobal::displayInfo(tdag.partialPathName() + tpoints[i].x + MString(" ") + tpoints[i].y + MString(" ") + tpoints[i].z + MString(" ") + tpoints[i].w);
		}
	}
	Z5Vector tempV = Z5Vector();
	tempV.Print();

	setResult("EmptyP command executed!\n");
	return MS::kSuccess;
}

MStatus YourselfCommand::undoIt()
{
	MGlobal::displayInfo("EmptyP command undone!\n");

	return MS::kSuccess;
}

bool YourselfCommand::isUndoable() const
{
	return true;//¿É³·ÏúµÄÃüÁî
}

void* YourselfCommand::creator()
{
	return new YourselfCommand();
}
