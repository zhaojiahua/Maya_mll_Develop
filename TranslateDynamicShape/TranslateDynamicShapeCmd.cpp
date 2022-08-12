#include "TranslateDynamicShapeCmd.h"

#include <maya/MGlobal.h>
#include <maya/MSelectionList.h>
#include <maya/MPoint.h>
#include <maya/MFnMesh.h>
#include <maya/MFloatArray.h>
#include <maya/MDoubleArray.h>
#include <maya/MArgList.h>


MStatus TranslateDynamicShape::doIt(const MArgList& arg)
{
	MStatus stat = MS::kSuccess;
	uvtolorance = arg.asDouble(0);
	return redoIt();
}

MStatus TranslateDynamicShape::redoIt()
{
	MSelectionList selections;
	MGlobal::getActiveSelectionList(selections);

	MStringArray selections_name;
	selections.getSelectionStrings(selections_name);
	MString baseName = selections_name[0] + "_Base";
	MStatus addState = selections.add(baseName);

	if (addState)
	{
		selections.getDagPath(0, dagPath1);
		selections.getDagPath(1, dagPath2);
		selections.getDagPath(2, dagPath3);

		MFnMesh fnMesh1(dagPath1);
		MFnMesh fnMesh2(dagPath2);
		MFnMesh fnMesh3(dagPath3);

		MStringArray mesh1uvnames, mesh2uvnames, mesh3uvnames;
		fnMesh1.getUVSetNames(mesh1uvnames);
		fnMesh2.getUVSetNames(mesh2uvnames);
		fnMesh3.getUVSetNames(mesh3uvnames);

		fnMesh2.getPoints(orgPoints);//首先获取mesh2的所有点当作原始点,可用于撤销命令

		MPointArray orgPoints_deff;//用来存储每个点的变形信息
		orgPoints_deff.clear();
		MIntArray pointMoved;//用来存储被移动的点
		
		for (unsigned i = 0; i < orgPoints.length();i++)
		{
			float2 tempuv;
			fnMesh2.getUVAtPoint(orgPoints[i], tempuv, MSpace::kObject, &mesh2uvnames[0]);

			int temppolygonID1, temppolygonID3;
			fnMesh1.intersectFaceAtUV(tempuv, temppolygonID1, &mesh1uvnames[0]);
			fnMesh3.intersectFaceAtUV(tempuv, temppolygonID3, &mesh3uvnames[0]);

			MPoint tempPoint1, tempPoint3;
			fnMesh1.getPointAtUV(temppolygonID1, tempPoint1, tempuv, MSpace::kWorld, &mesh1uvnames[0], uvtolorance);
			fnMesh3.getPointAtUV(temppolygonID3, tempPoint3, tempuv, MSpace::kWorld, &mesh3uvnames[0], uvtolorance);

			MVector PointDeff = tempPoint1 - tempPoint3;
			
			if (PointDeff.length()>0.001)
			{
				pointMoved.append(i);
			}
			orgPoints_deff.append(orgPoints[i] + PointDeff);
		}
		fnMesh2.setPoints(orgPoints_deff);

		setResult(pointMoved);
	}
	else 
	{
		displayError("There is no mesh named " + baseName + "!");
	}

	return MS::kSuccess;
}

MStatus TranslateDynamicShape::undoIt()
{
    MGlobal::displayInfo( "TranslateDynamicShape command undone!\n" );

	MFnMesh fnMesh2(dagPath2);
	fnMesh2.setPoints(orgPoints, MSpace::kObject);

	return MS::kSuccess;
}

void* TranslateDynamicShape::creator()
{
	return new TranslateDynamicShape();
}

TranslateDynamicShape::TranslateDynamicShape()
{}

TranslateDynamicShape::~TranslateDynamicShape()
{}

bool TranslateDynamicShape::isUndoable() const
{
	return true;
}
