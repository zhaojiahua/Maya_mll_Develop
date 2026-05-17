#include "ImportSkinClusterDatasCmd.h"
#include <maya/MGlobal.h>

ImportSkinClusterDatas::ImportSkinClusterDatas():m_filePathStr("C:/mayaSkinWeights/") {}

MStatus ImportSkinClusterDatas::getSkinCluster(MObject inNode, MObject& skinCluster){
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

MSyntax ImportSkinClusterDatas::newSyntax() {
	MSyntax syntax;
	syntax.addFlag("-p", "-path", MSyntax::kString);
	return syntax;
}

MStatus ImportSkinClusterDatas::parseArgs(const MArgList& args){
	MArgDatabase argData(syntax(), args);
	if (argData.isFlagSet("-p")) {
		argData.getFlagArgument("-p", 0, m_filePathStr);
	}
	else {
		displayError("----->>请指定.csv的文件夹路径!");
		return MS::kFailure;
	}
	return MS::kSuccess;
}

MStatus ImportSkinClusterDatas::readCSVFile(const MString& fullpath){
	ifstream file(fullpath.asChar());
	if (!file.is_open()) {
		displayError(MString("无法打开文件:") + fullpath);
		return MS::kFailure;
	}
	string line;
	//读取首行
	{
		if (!std::getline(file, line)) {
			displayError(".csv文件为空!");
			return MS::kFailure;
		}
		//存储读取到的joints名称
		m_jointsNames.clear();
		stringstream strstream(line);
		string token;
		getline(strstream, token, ',');//跳过首个元素-顶点序号
		while (getline(strstream, token, ',')) {
			m_jointsNames.append(token.c_str());
		}
	}
	//读取顶点权重
	{
		while (getline(file,line)){
			stringstream lineStrStream(line);
			string weightStr;
			getline(lineStrStream, weightStr, ',');//跳过首个元素
			MDoubleArray tempWeights;
			for (auto& jointname : m_jointsNames) {
				getline(lineStrStream, weightStr, ',');
				tempWeights.append(stod(weightStr));
			}
			m_vertexWeights.push_back(tempWeights);
		}
	}
	file.close();
	return MStatus();
}

MStatus ImportSkinClusterDatas::createSkinCluster(){
	if (m_skinClusterNode_old != MObject::kNullObj) {
		m_skinClusterNode = m_skinClusterNode_old;
		return MS::kSuccess;
	}
	MDagPathArray inflenceJointsPaths;
	for (const auto& jname : m_jointsNames) {
		MSelectionList tempSl;
		tempSl.add(jname);
		MDagPath tempJoinPath;
		if (tempSl.getDagPath(0, tempJoinPath) == MS::kSuccess) {
			inflenceJointsPaths.append(tempJoinPath);
		}
		else {
			displayWarning(jname + "在场景中没有找到!");
		}
	}
	if (inflenceJointsPaths.length() == 0) {
		displayError("没有找到有效的joint");
		return MS::kFailure;
	}
	MFnSkinCluster skinFn{};
	MDGModifier dgModifier;
	MDoubleArray initVertexWeights{};
	for (auto joint : inflenceJointsPaths)initVertexWeights.append(1.0);
	MObject skinNode = skinFn.create(m_meshDagPath, inflenceJointsPaths, initVertexWeights, dgModifier);

	return MStatus();
}

MStatus ImportSkinClusterDatas::applyWeights(MObject& inSkinCluster, const vector<MDoubleArray>& inVertexWeights, vector<MDoubleArray>& inVertexWeights_old) {
	if (inSkinCluster == MObject::kNullObj || inSkinCluster.apiType() != MFn::kSkinClusterFilter)return MS::kFailure;
	MFnSkinCluster skinFn(inSkinCluster);
	//获取权重列表中这些joints的排序
	MIntArray influenceIndices;
	for (unsigned int i = 0; i < m_jointsNames.length(); ++i) {
		MStatus status;
		MSelectionList tempsl{};
		tempsl.add(m_jointsNames[i]);
		MDagPath tempJointDagPath{};
		tempsl.getDagPath(0, tempJointDagPath);
		unsigned int idx = skinFn.indexForInfluenceObject(tempJointDagPath, &status);	//获取joint的排序序号
		if (status == MS::kSuccess) {
			influenceIndices.append((int)idx);
		}
		else {
			// Fallback: assume indices are in order 0..N-1
			influenceIndices.append(i);
		}
	}
	for (unsigned int v = 0; v < inVertexWeights.size(); ++v) {
		MDoubleArray vertexWeights{};	//按照joints的顺序,构建权重数据
		for (unsigned int i = 0; i < influenceIndices.length(); ++i) {
			vertexWeights.append(inVertexWeights.at(v)[i]);
		}
	 	const MStatus stat = skinFn.setWeights(m_meshDagPath, m_meshComponent, influenceIndices, vertexWeights, false, &(inVertexWeights_old[v]));
		if (stat == MS::kSuccess) {
			displayWarning(MString("----->>顶点") + v + "的权重设置失败!");
		}
	}
	displayInfo(MString("-------------------------------------->>") + m_meshDagPath.partialPathName() + "蒙皮权重设置成功!");
	return MStatus();
}

//先选目标体再选本体
MStatus ImportSkinClusterDatas::doIt( const MArgList& arglist){
	const MStatus stat = parseArgs(arglist);
	if (stat != MS::kSuccess) {
		displayError("----->>参数解析失败!");
		return stat; 
	}
	return redoIt();
}

MStatus ImportSkinClusterDatas::redoIt()
{
	MGlobal::getActiveSelectionList(m_selectionlist);
	if (m_selectionlist.length() != 1) {
		displayError("----->>必须选择且只能选择一个mesh物体");
		return MS::kFailure;
	}
	MStatus stat = m_selectionlist.getDagPath(0, m_meshDagPath);
	if (stat != MS::kSuccess) {
		displayError("----->>获取所选物体的DAG Path失败!");
		return stat;
	}
	stat = readCSVFile(m_filePathStr + m_meshDagPath.partialPathName() + ".csv");
	if (stat != MS::kSuccess) {
		displayError("----->>readCSVFile执行失败!");
		return stat;
	}
	//先查找其原有的skinCluster节点
	stat = getSkinCluster(m_meshDagPath.node(), m_skinClusterNode_old);
	if (stat != MS::kSuccess)m_skinClusterNode_old = MObject::kNullObj;
	//再创建新的skinCluster节点
	stat = createSkinCluster();
	if (stat != MS::kSuccess) {
		displayError("----->>createSkinCluster执行失败!");
		return stat;
	}
	stat = applyWeights(m_skinClusterNode, m_vertexWeights, m_vertexWeights_old);
	if (stat != MS::kSuccess) {
		displayError("----->>applyWeights应用权重数据失败!");
		return stat;
	}
	return MS::kSuccess;
}

MStatus ImportSkinClusterDatas::undoIt(){
	if (m_skinClusterNode_old == MObject::kNullObj) {
		MDGModifier dagModifier;
		dagModifier.deleteNode(m_skinClusterNode);
		return dagModifier.doIt();
	}
	else {
		const MStatus stat = applyWeights(m_skinClusterNode_old, m_vertexWeights_old, m_vertexWeights);
		if (stat != MS::kSuccess) {
			displayWarning("权重恢复失败!");
			return MS::kFailure;
		}
		displayInfo("撤销操作成功,权重已恢复为原来的数据");
	}
	return MS::kSuccess;
}

void* ImportSkinClusterDatas::creator() { return new ImportSkinClusterDatas(); }

ImportSkinClusterDatas::~ImportSkinClusterDatas()
{}

bool ImportSkinClusterDatas::isUndoable() const
{
	return true;//可撤销命令
}
