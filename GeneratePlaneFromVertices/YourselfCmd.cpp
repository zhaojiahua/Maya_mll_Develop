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
	for (unsigned int i = 0; i < allselections.length(); ++i) {
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
		suc = tfnmesh.getPoints(tpoints, MSpace::kWorld);

	}
	Z5Matrix newmat = Z5Matrix();
	for (int i = 0; i < 5; ++i) {
		for (int j = 0; j < 5; ++j) {
			newmat.SetElement(i, j, (int)(MRandom::Rand_d((i + 1) * (j + 1) + 17, 113) * 100));
		}
	}
	newmat.Print();
	Z5Matrix inversemat = newmat.InverseMatrix();
	inversemat.Print();
	Z5Matrix nresult = newmat * inversemat;
	nresult.Print();

	return MS::kSuccess;
}

MStatus YourselfCommand::undoIt()
{
	MGlobal::displayInfo("EmptyP command undone!\n");

	return MS::kSuccess;
}

bool YourselfCommand::isUndoable() const
{
	return true;
}

void* YourselfCommand::creator()
{
	return new YourselfCommand();
}

MStatus YourselfCommand::Get2DPoints(const MDagPath& inpath, int&, MDoubleArray&, double&)
{
	MStatus succ;
	MFnMesh tfnmesh(inpath);
	MPointArray allpoints;
	succ = tfnmesh.getPoints(allpoints, MSpace::kWorld);
	unsigned int ptnum = tfnmesh.numVertices();
	MDoubleArray xdatas ;
	MDoubleArray ydatas ;
	MDoubleArray zdatas ;
	MDoubleArray xzdatas ;//x轴
	MDoubleArray xydatas ;//x轴
	MDoubleArray yxdatas ;//y轴
	MDoubleArray yzdatas ;//y轴
	MDoubleArray zxdatas ;//z轴
	MDoubleArray zydatas ;//z轴
	for (int i = 0; i < ptnum; ++i) {
		xdatas.append(allpoints[i].x);
		ydatas.append(allpoints[i].y);
		zdatas.append(allpoints[i].z);
		xzdatas.append(allpoints[i].x); xzdatas.append(allpoints[i].z);
		xydatas.append(allpoints[i].x); xydatas.append(allpoints[i].y);
		yxdatas.append(allpoints[i].y); xydatas.append(allpoints[i].x);
		yzdatas.append(allpoints[i].y); xydatas.append(allpoints[i].z);
		zxdatas.append(allpoints[i].z); xydatas.append(allpoints[i].x);
		zydatas.append(allpoints[i].z); xydatas.append(allpoints[i].y);
	}
	
	return succ;
}
