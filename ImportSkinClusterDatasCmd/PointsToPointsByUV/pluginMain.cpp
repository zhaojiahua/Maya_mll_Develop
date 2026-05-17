#include "ImportSkinClusterDatasCmd.h"
#include <maya/MFnPlugin.h>

MStatus initializePlugin( MObject obj )
{ 
	MStatus   status;
	MFnPlugin plugin( obj, "ZhaoJiahua", "2022", "Any");

	status = plugin.registerCommand( "ImportSkinClusterDatas", ImportSkinClusterDatas::creator );
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

	status = plugin.deregisterCommand( "ImportSkinClusterDatas" );
	if (!status) {
		status.perror("deregisterCommand");
		return status;
	}

	return status;
}
