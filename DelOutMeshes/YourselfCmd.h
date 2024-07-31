#pragma once
#include <maya/MPxCommand.h>
#include <maya/MSelectionList.h>
#include <maya/MFnTransform.h>
#include <maya/MDagPath.h>
#include <maya/MDagPathArray.h>
#include <maya/MPoint.h>
#include <maya/MPointArray.h>
#include <maya/MFnMesh.h>
#include <maya/MVector.h>
#include <maya/MIntArray.h>
#include<maya/MArgList.h>
#include <maya/MArgParser.h>
#include <maya/MSyntax.h>
#include <maya/MArgParser.h>


//Mayaµƒ√¸¡Ó¿‡
class YourselfCommand :public MPxCommand
{
public:
	YourselfCommand();
	virtual		~YourselfCommand();

	MStatus		doIt(const MArgList& args);
	MStatus		redoIt();
	MStatus		undoIt();
	bool		isUndoable() const;

	static		void* creator();

private:
	// Store the data you will need to undo the command here
	//
	MSelectionList sele_ls;
	MDagPathArray sele_dags;
	bool issensitive = false;
};