#include "ExportSkinClusterDatasCmd.h"
#include <maya/MGlobal.h>

ExportSkinClusterDatas::ExportSkinClusterDatas():filePathStr("C:/mayaSkinWeights/") {}

MStatus ExportSkinClusterDatas::getSkinCluster(const MObject& inNode, MObject& skinCluster){
	MFnDagNode dagNode(inNode);
	const unsigned int childcount= dagNode.childCount();
	for (unsigned int i = 0; i < childcount; ++i) {
		MObject childMesh = dagNode.child(i);
		if(childMesh.apiType() == MFn::kMesh) {
			MFnDagNode kMeshNode(childMesh);
			MStatus findInMeshPlugStat;
			MPlug inMeshPlug = kMeshNode.findPlug("inMesh", findInMeshPlugStat);
			if (findInMeshPlugStat == MS::kSuccess) {
				MObject inMeshSourceNode = inMeshPlug.source().node();
				MFnDagNode inMeshSourceDagNode(inMeshSourceNode);
				MString tempForOutStr = kMeshNode.name() + "的inMesh接口有连接,它的连接源是" + inMeshSourceDagNode.name();
				if (inMeshSourceNode.apiType() == MFn::kSkinClusterFilter) {
					skinCluster = inMeshSourceNode;
					displayInfo("----->> 通过inMesh接口查找到了SkinClusterFilter节点" + inMeshSourceDagNode.name());
					return MS::kSuccess;
				}
				else 	displayWarning(kMeshNode.name() + "的inMesh接口的连接源头不是kSkinClusterFilter节点!");
			}
			else displayWarning("---->> 未找到" + kMeshNode.name() + "的inMesh接口!");
			//如果通过inMesh接口没有查找到SkinCluster节点,就查找其skinClusterSet节点,然后通过skinClusterSet获取skinCluster节点
			MPlugArray plugs;
			kMeshNode.getConnections(plugs);
			for (auto& plug : plugs) {
				if (plug.isDestination()) {
					MObject sourceNode = plug.source().node();
					if (sourceNode.apiType() == MFn::kSet) {
						MFnDependencyNode sourceDependencyNode(sourceNode);
						MPlugArray setPlugs;
						sourceDependencyNode.getConnections(setPlugs);
						for (auto& setPlug : setPlugs) {
							if (setPlug.isDestination()) {
								MObject setSourceNode = setPlug.source().node();
								if (setSourceNode.apiType() == MFn::kSkinClusterFilter) {
									skinCluster = setSourceNode;
									MFnDependencyNode finalSkinNode(setSourceNode);
									displayInfo("----->> 通过skinClusterSet查找到了SkinClusterFilter节点" + finalSkinNode.name());
									return MS::kSuccess;
								}
							}
						}
					}
				}
			}
		}
	}
	return MS::kFailure;
}

MStatus ExportSkinClusterDatas::getSkinCluster_new(MObject inNode, MObject& skinCluster){
	//使用MItDependencyGraph进行查找节点的方法,相比于上一个方法,这个方法不需要知道具体的连接关系,只需要指定查找的节点类型和连接方向就可以了
	MFnDagNode dagNode(inNode);
	const unsigned int childcount = dagNode.childCount();
	for (unsigned int i = 0; i < childcount; ++i) {
		MObject childMesh = dagNode.child(i);
		if (childMesh.apiType() == MFn::kMesh) {
			MItDependencyGraph dgIt_upper(childMesh, MFn::kSkinClusterFilter, MItDependencyGraph::kUpstream);
			if (!dgIt_upper.isDone()) {
				skinCluster = dgIt_upper.currentItem();
				return MS::kSuccess;
			}
		}
	}
	return MS::kFailure;
}

MStatus ExportSkinClusterDatas::parseArgs(const MArgList& args){
	//强制输入参数 flag: -f/-file <filename>
	MStatus stat;
	MString arg;
	const MString fileFlag("-fp");
	const MString fileFlagLong("-filePath");
	//解析参数
	for (unsigned int i = 0; i < args.length(); ++i) {
		arg = args.asString(i, &stat);
		if (!stat) continue;
		if(arg == fileFlag || arg == fileFlagLong){
			//获取传进来的文件路径参数
			if(i==args.length()-1){
				arg += "---->>必须指定文件夹路径!";
				displayError(arg);	return MS::kFailure;
			}
			i++;
			args.get(i, filePathStr);
		}
		else {
			arg += "---->> 未知参数!";
			displayError(arg); return MS::kFailure;
		}
	}
	return MS::kSuccess;
}

//先选目标体再选本体
MStatus ExportSkinClusterDatas::doIt( const MArgList& arglist)
{
	//解析参数,获取文件路径
	MStatus stat = parseArgs(arglist);
	if (stat != MS::kSuccess) { return stat; }

	unsigned int skincount = 0;
	//遍历所选择的节点,找到所有的skinCluster节点
	MSelectionList selList;
	MGlobal::getActiveSelectionList(selList);
	MItSelectionList selIter(selList, MFn::kInvalid);
	for (; !selIter.isDone(); selIter.next()) {
		MObject selectedDependNode;
		selIter.getDependNode(selectedDependNode);
		MFnDependencyNode fnDependencyNode(selectedDependNode);
		MObject skinClusterNode;//skinCluster节点
		stat = getSkinCluster_new(selectedDependNode, skinClusterNode);
		if (stat == MS::kSuccess) {
			skincount++;
			//对于每个skinCluster节点,先获取它的所有影响物体
			MFnSkinCluster skinCluster(skinClusterNode);
			//所有的蒙皮权重影响体
			MDagPathArray infObjPathes;
			const unsigned int nInfs = skinCluster.influenceObjects(infObjPathes, &stat);
			CheckError(stat, "---->>获取影响物体失败!");
			if (nInfs == 0) {
				stat = MS::kFailure;
				CheckError(stat, "---->>没有权重物体");
			}

			//打开文件,准备写入
			MString fileFullPathStr = filePathStr + fnDependencyNode.name() + "_skinWeights.csv";
			FILE* file = fopen(fileFullPathStr.asChar(), "wb");
			if (!file) {
				CheckError(MS::kFailure, "---->>无法打开文件:" + fileFullPathStr);
			}
			//遍历这个skinCluster节点的所有影响Mesh
			const unsigned int nGeoms = skinCluster.numOutputConnections();
			for (unsigned int j = 0; j< nGeoms; j++) {
				const unsigned int index = skinCluster.indexForOutputConnection(j, &stat);
				CheckError(stat, "---->>获取Geometry index 失败!");
				MDagPath skinPath;
				stat = skinCluster.getPathAtIndex(index, skinPath);
				CheckError(stat, "---->>获取Geometry path 失败!");
				//遍历Geometry的所有Components
				MItGeometry GeometryIter(skinPath);//可以用来遍历polygon的vertices

				//打印出所有影响体
				fprintf(file, " ,");//开头第一个空格
				for (unsigned int i= 0; i < nInfs; ++i) {
					fprintf(file, "%s,", infObjPathes[i].partialPathName().asChar());
				}
				fprintf(file, "\n");
				unsigned int progressCount{ 0 };
				const unsigned int processStep = 0.05f * GeometryIter.count();//每5%打印一次
				int progress{ 0 };
				for (; !GeometryIter.isDone(); GeometryIter.next()) {
					MObject comp = GeometryIter.currentItem(&stat);
					CheckError(stat, "---->>获取Component失败!");
					//获取每个点的权重信息
					MDoubleArray weights;
					unsigned int infCount;
					stat = skinCluster.getWeights(skinPath, comp, weights, infCount);
					CheckError(stat, "---->>获取权重失败!");
					if (infCount == 0) {
						stat = MS::kFailure;
						CheckError(stat, "---->>没有权重!");
					}
					//输出这些点的权重数据
					fprintf(file, "%d,", GeometryIter.index());
					for (auto& weight : weights) {
						fprintf(file, "%f,", weight);
					}
					fprintf(file, "\n");
					progressCount++;
					if (progressCount > processStep) {
						progress += 5;
						displayInfo(MString("==========>>") + progress + "%");
						progressCount = 0;
					}
				}
			}
			displayInfo("------------------------------=============>>100%");
			displayInfo(fnDependencyNode.name() + "权重数据成功写入文件:" + fileFullPathStr);
			fclose(file);
		}
		else { CheckError(MS::kFailure, fnDependencyNode.name() + "获取skinCluster节点失败!"); }
	}
	return MS::kSuccess;
}

MStatus ExportSkinClusterDatas::redoIt()
{
	clearResult();
	setResult(1);
	return MS::kSuccess;
}

MStatus ExportSkinClusterDatas::undoIt()
{
    MGlobal::displayInfo( "PointsToPointsByUV command undone!\n" );

	return MS::kSuccess;
}

void* ExportSkinClusterDatas::creator()
{
	return new ExportSkinClusterDatas();
}


ExportSkinClusterDatas::~ExportSkinClusterDatas()
{}

bool ExportSkinClusterDatas::isUndoable() const
{
	return false;//不可撤销命令
}
