#ifndef _TranslateDynamicShapeCmd
#define _TranslateDynamicShapeCmd

#include <maya/MPxCommand.h>
#include <maya/MDagPath.h>
#include <maya/MPointArray.h>

class MArgList;

class TranslateDynamicShape : public MPxCommand
{

public:
				TranslateDynamicShape();
	virtual		~TranslateDynamicShape();

	MStatus		doIt( const MArgList& arg );
	MStatus		redoIt();
	MStatus		undoIt();
	bool		isUndoable() const;

	static		void* creator();

private:
	// Store the data you will need to undo the command here
	MPointArray orgPoints;

	MDagPath dagPath1, dagPath2, dagPath3;

	double uvtolorance = 0.01;
};

#endif
