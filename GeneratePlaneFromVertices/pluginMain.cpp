#include "YourselfCmd.h"
#include <maya/MFnPlugin.h>

MStatus initializePlugin(MObject obj)
{
	MStatus status;
	MFnPlugin plugin(obj, "Zjh", "2022", "Any");//指明该插件的作者 版本 和使用的API
	status = plugin.registerCommand("GeneratePlaneFromVertices", YourselfCommand::creator);//这里向Maya注册命令,提供启用该命令的字符串名称
	if (!status) {
		status.perror("registerCommand");
		return status;
	}

	return status;
}
MStatus uninitializePlugin(MObject obj)
{
	MStatus   status;
	MFnPlugin plugin(obj);

	status = plugin.deregisterCommand("GeneratePlaneFromVertices");
	if (!status) {
		status.perror("deregisterCommand");
		return status;
	}

	return status;
}