#include "ExportSkinClusterDatasCmd.h"
#include <maya/MFnPlugin.h>

MStatus initializePlugin( MObject obj )
{ 
	MStatus   status;
	MFnPlugin plugin( obj, "ZhaoJiahua", "2020", "Any");

	status = plugin.registerCommand( "ExportSKinClusterDatas", ExportSkinClusterDatas::creator );
	if (!status) {
		status.perror("registerCommand");
		return status;
	}

	return status;
}

MStatus uninitializePlugin( MObject obj )
{
	MStatus   status;
	MFnPlugin plugin( obj );

	status = plugin.deregisterCommand( "ExportSKinClusterDatas" );
	if (!status) {
		status.perror("deregisterCommand");
		return status;
	}

	return status;
}
