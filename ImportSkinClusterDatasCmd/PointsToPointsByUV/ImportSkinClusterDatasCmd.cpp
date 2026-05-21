#include "ImportSkinClusterDatasCmd.h"
#include <maya/MGlobal.h>

ImportSkinClusterDatas::ImportSkinClusterDatas():filePath("C:/mayaSkinWeights/") {}

MSyntax ImportSkinClusterDatas::newSyntax() {
	MSyntax syntax;
	syntax.addFlag("-fp", "-filePath", MSyntax::kString);
	return syntax;
}

MStatus ImportSkinClusterDatas::parseArgs(const MArgList& args){
	MArgDatabase argData(newSyntax(), args);
	if (argData.isFlagSet("-fp")) {
		argData.getFlagArgument("-fp", 0, filePath);
	}
	else {
		displayError("----->>请指定.csv的文件夹路径!");
		return MS::kFailure;
	}
	return MS::kSuccess;
}

MStatus ImportSkinClusterDatas::GetSkinCluster(MDagPath inPath, MObject& skinCluster){
    if (inPath.node().apiType() != MFn::kMesh) {
        MStatus stat = inPath.extendToShape();
        if (stat != MS::kSuccess) {
            displayError(inPath.partialPathName() + "---->> 不是Mesh类型,并且无法找到它的Shape节点!");
            return MS::kFailure;
        }
    }
    MItDependencyGraph dgIt_upper(inPath.node(), MFn::kSkinClusterFilter, MItDependencyGraph::kUpstream);
    if (!dgIt_upper.isDone()) {
        skinCluster = dgIt_upper.currentItem();
        return MS::kSuccess;
    }
    return MS::kFailure;
}

MStatus ImportSkinClusterDatas::ReadCSVFile(const MString& fullpath, MStringArray& m_jointsNames, vector<MDoubleArray>& m_vertexWeights){
    std::ifstream file(fullpath.asChar());
    if (!file.is_open()) {
        MString msg = "Cannot open CSV file: ";
        MGlobal::displayError(msg + fullpath);
        return MS::kFailure;
    }
    std::string line;
    // Read header line: vertex_index,joint1,joint2,joint3,...
    if (!std::getline(file, line)) {
        MGlobal::displayError(fullpath + " is empty");
        return MS::kFailure;
    }
    // Parse joint names
    std::stringstream ss(line);
    std::string tempjoint;
    std::getline(ss, tempjoint, ','); // skip "vertex_index"
    while (std::getline(ss, tempjoint, ',')) {
        m_jointsNames.append(MString(tempjoint.c_str()));
    }
    // Read weight data
    while (std::getline(file, line)) {
        std::stringstream lineStream(line);
        std::string tempweight;
        // skip vertex index
        std::getline(lineStream, tempweight, ',');
        MDoubleArray weights;
        for (unsigned int i = 0; i < m_jointsNames.length(); ++i) {
            std::getline(lineStream, tempweight, ',');
            weights.append(std::stod(tempweight));
        }
        m_vertexWeights.push_back(weights);
    }
    file.close();
    return MS::kSuccess;
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
    MStatus stat;
    MSelectionList templist{};
    MGlobal::getActiveSelectionList(templist);
    for (unsigned int i = 0; i < templist.length(); ++i) {
        MDagPath tempDagPath;
        templist.getDagPath(i, tempDagPath);
        //if (tempDagPath.apiType() != MFn::kMesh)tempDagPath.extendToShape();
        selectedDagPaths.append(tempDagPath);
    }
    MProgressWindow::reserve();
    MProgressWindow::setTitle("正在导入权重数据!");
    MProgressWindow::setInterruptable(true);
    MProgressWindow::startProgress();
    unsigned int processIndex = 0;
    for (auto selectedDagPath : selectedDagPaths) {
        MStringArray m_jointsNames{};//存储骨骼名称
        vector<MDoubleArray> m_vertexWeights{};//存储定点权重
        stat = ReadCSVFile(filePath + selectedDagPath.partialPathName() + "_skinWeights.csv", m_jointsNames, m_vertexWeights);
        MItGeometry GeoIt(selectedDagPath);
        MIntArray vertexIndices{};
        for (; !GeoIt.isDone(); GeoIt.next()) {
            vertexIndices.append(GeoIt.index());
        }
        if (m_vertexWeights.size() != vertexIndices.length()) {
            displayError(selectedDagPath.partialPathName() + "读取的.csv权重数目与其vertexComponent数目不一致! 读取的行数:" + m_vertexWeights.size() + "自身的Mesh顶点数目:" + vertexIndices.length());
            continue;
        }
        //创建skinCluster
        if (stat == MS::kSuccess) {
            MString melCmdLine("skinCluster -dr 10 -tsb ");
            for (auto jointname : m_jointsNames) {
                melCmdLine += jointname + " ";
            }
            melCmdLine += selectedDagPath.partialPathName();
            dgModifier.commandToExecute(melCmdLine);
            dgModifier.doIt();
            //获取skinCluster节点
            MObject skinClusterNode;
            stat = GetSkinCluster(selectedDagPath, skinClusterNode);
            //应用权重
            if (stat == MS::kSuccess) {
                MFnSkinCluster skinFn(skinClusterNode);
                //构建骨骼顺序参数
                MIntArray inflenceIndices{};
                for (auto jointname : m_jointsNames) {
                    MSelectionList tempsl{};
                    tempsl.add(jointname);
                    MDagPath joinDagPath{};
                    tempsl.getDagPath(0, joinDagPath);
                    const unsigned int jointInflenceIndex = skinFn.indexForInfluenceObject(joinDagPath);
                    inflenceIndices.append(jointInflenceIndex);
                }
                //构建权重参数
                MDoubleArray tempweights;
                for (auto vertxweight : m_vertexWeights) {
                    for (auto dweight : vertxweight) {
                        tempweights.append(dweight);
                    }
                }
                //构建component参数
                MFnSingleIndexedComponent compFn;
                MObject allverticesComp = compFn.create(MFn::kMeshVertComponent);
                compFn.addElements(vertexIndices);
                //应用权重
                stat = skinFn.setWeights(selectedDagPath, allverticesComp, inflenceIndices, tempweights, false);
                if (stat == MS::kSuccess) {
                    displayInfo(selectedDagPath.partialPathName() + "--------------------------------------->>权重应用成功!");
                }
            }
            else displayError(selectedDagPath.partialPathName() + "GetSkinCluster失败!");
        }
        processIndex++;
        MProgressWindow::setProgress(static_cast<float>(processIndex) / static_cast<float>(selectedDagPaths.length()) * 100);
    }
    MProgressWindow::endProgress();
    return MS::kSuccess;
}

MStatus ImportSkinClusterDatas::undoIt() {
    dgModifier.undoIt();
    displayInfo("ImportSkinClusterDatas Command undone!\n");
    return MS::kSuccess;
}

void* ImportSkinClusterDatas::creator() { return new ImportSkinClusterDatas(); }

ImportSkinClusterDatas::~ImportSkinClusterDatas()
{}

bool ImportSkinClusterDatas::isUndoable() const
{
	return true;//可撤销命令
}
