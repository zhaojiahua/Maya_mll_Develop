#include <maya/MSimple.h>
#include <maya/MFnMesh.h>
#include <maya/MIntArray.h>
#include <maya/MDoubleArray.h>
#include <maya/MPointArray.h>
#include <maya/MPoint.h>
#include <maya/MStringArray.h>
#include <maya/MString.h>
#include <maya/MSelectionList.h>
#include <maya/MGlobal.h>
#include <maya/MDagPath.h>

DeclareSimpleCommand( GetPointsAtUV, "zhaoJiaHua", "2018");

//根据输入的UV值,获取该面上的距离该uv点最近的面和点,返回faceID和vertexID

MStatus GetPointsAtUV::doIt( const MArgList& args )
{
	MStatus stat = MS::kSuccess;
	MSelectionList selection;
	stat = MGlobal::getActiveSelectionList(selection);
	if (selection.length()==1)
	{
		MDagPath path1;
		selection.getDagPath(0, path1);
		MFnMesh fnMesh1(path1);
		MIntArray polygonIDs;
		for (int i = 0; i < fnMesh1.numPolygons(); i++)
		{
			polygonIDs.append(i);
		}

		MStringArray uvnames;
		fnMesh1.getUVSetNames(uvnames);

		double uvalue = args.asDouble(0);
		double vvalue = args.asDouble(1);
		float2 uvvalue = { uvalue, vvalue };
		MPointArray points;

		int outPolygonID = 0;
		fnMesh1.intersectFaceAtUV(uvvalue, outPolygonID, &uvnames[0]);//给定UVset和UV坐标,返回这个VU坐标下所在的FaceID
		MPoint uvPoint;
		fnMesh1.getPointAtUV(outPolygonID, uvPoint, uvvalue, MSpace::kWorld, &uvnames[0]);//获得指定UVset下指定UV值的世界空间位置


		MIntArray polygonVertices;
		fnMesh1.getPolygonVertices(outPolygonID, polygonVertices);

		MDoubleArray distanceArray;
		distanceArray.clear();
		for (unsigned i = 0; i < polygonVertices.length();i++)
		{
			MPoint tempPoint;
			fnMesh1.getPoint(polygonVertices[i], tempPoint, MSpace::kWorld);
			distanceArray.append(tempPoint.distanceTo(uvPoint));
		}

		//冒泡排序查找最近的vertex
		for (unsigned i = 0; i < distanceArray.length() - 1; i++)
		{
			for (unsigned j = 0; j < distanceArray.length() - 1 - i; j++)
			{
				if (distanceArray[j]>distanceArray[j+1])
				{
					double tempdis = distanceArray[j];
					distanceArray[j] = distanceArray[j + 1];
					distanceArray[j + 1] = tempdis;
					int tempid = polygonVertices[j];
					polygonVertices[j] = polygonVertices[j + 1];
					polygonVertices[j + 1] = tempid;
				}
			}
		}

		MDoubleArray outResultArray;
		outResultArray.clear();
		outResultArray.append(outPolygonID);//查找到的FaceID传出去
		outResultArray.append(polygonVertices[0]);//查找到的最近的VertexID传出去
		//把查找到的世界空间位置传出去
		outResultArray.append(uvPoint.x);
		outResultArray.append(uvPoint.y);
		outResultArray.append(uvPoint.z);
		setResult(outResultArray);
	}
	else { MGlobal::displayError("this is no mesh was selected!"); }
	return stat;
}
