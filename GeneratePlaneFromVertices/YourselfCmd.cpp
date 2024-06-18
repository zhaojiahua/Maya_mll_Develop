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
	for (unsigned int i=0; i< allselections.length();++i) {
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
		//for (int i = 0; i < tpoints.length(); ++i){
		//	MGlobal::displayInfo(tdag.partialPathName() + tpoints[i].x + MString(" ") + tpoints[i].y + MString(" ") + tpoints[i].z + MString(" ") + tpoints[i].w);
		//}
	}
	Z5Matrix newmat = Z5Matrix();
	for (int i = 0; i < 5; ++i) {
		for (int j = 0; j < 5; ++j) {
			newmat.SetElement(i, j, (int)(MRandom::Rand_d((i + 1) * (j + 1) + 13, 11) * 100));
		}
	}
	newmat.Print();
	MDoubleArray submat1 = Z5Matrix::SubMatrix(newmat.GetMateData(), 1, 2);
	MDoubleArray submat2 = Z5Matrix::SubMatrix(submat1, 1,1);
	Z5Matrix::PrintMateDatas(submat1);
	Z5Matrix::PrintMateDatas(submat2);


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
