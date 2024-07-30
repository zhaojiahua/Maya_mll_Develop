#include "YourselfCmd.h"
#include <utility>
#include <maya/MGlobal.h>
#include <maya/MSyntax.h>
#include <maya/MTime.h>
#include <maya/MAnimControl.h>
#include <maya/MArgParser.h>

YourselfCommand::YourselfCommand()
{}
YourselfCommand::~YourselfCommand()
{}

MStatus YourselfCommand::doIt(const MArgList& args)
{

	MGlobal::getActiveSelectionList(sele_ls);
	for (int i = 0; i < sele_ls.length(); ++i) {
		MDagPath tempdag;
		sele_ls.getDagPath(i, tempdag);
		sele_dags.append(tempdag);
	}
	return redoIt();
}

MStatus YourselfCommand::redoIt()
{
	//MFnTransform fnTrans(sele_dags[1]);
	//MPoint tworld = fnTrans.rotatePivot(MSpace::kWorld);
	MFnMesh fnMehs0(sele_dags[0]);
	MFnMesh fnMehs1(sele_dags[1]);
	int faceCount = fnMehs1.numPolygons();
	MPointArray allpoints;
	fnMehs1.getPoints(allpoints, MSpace::kWorld);
	MIntArray facesDelete;
	for (int i = 0; i < faceCount; ++i) {
		MIntArray vertexlist;
		fnMehs1.getPolygonVertices(i, vertexlist);
		int j = 0;
		for (; j < vertexlist.length(); ++j) {
			MPoint selfPoint = allpoints[vertexlist[j]];
			MPoint closetPosition;
			MVector closetNormal;
			int closetFaceId;
			fnMehs0.getClosestPointAndNormal(selfPoint, closetPosition, closetNormal, MSpace::kWorld, &closetFaceId, NULL);
			MVector direc = selfPoint - closetPosition;
			if (direc * closetNormal < 0) break;	
		}
		if (j == vertexlist.length())facesDelete.append(i);
	}
	if (facesDelete.length() == faceCount) {
		facesDelete = MIntArray();
		facesDelete.append(0);
	}
	setResult(facesDelete);
	//MGlobal::displayInfo(MString("closetPt") + closetPt.x + "," + closetPt.y + "," + closetPt.z);
	//MGlobal::displayInfo(MString("polynormal") + polynormal.x + "," + polynormal.y + "," + polynormal.z);
	//MGlobal::displayInfo(MString("polynum: ") + fnMehs.numPolygons());
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
