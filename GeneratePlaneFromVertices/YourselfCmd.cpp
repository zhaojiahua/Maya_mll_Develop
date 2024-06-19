#include "YourselfCmd.h"


YourselfCommand::YourselfCommand()
{}
YourselfCommand::~YourselfCommand()
{}

MStatus YourselfCommand::doIt(const MArgList& arglist)
{
	MStatus stat = MS::kSuccess;
	MSelectionList allselections;
	MGlobal::getActiveSelectionList(allselections);
	MStatus suc;
	for (unsigned int i = 0; i < allselections.length(); ++i) {
		MDagPath tdag;
		allselections.getDagPath(i, tdag);
		MFnDagNode tfnnode(tdag);
		if (tfnnode.child(0).apiType() == MFn::kMesh) {
			selections.add(tdag);
		}
	}
	MSyntax syt = MSyntax();
	syt.addFlag("cs", "curveSpans", MSyntax::kUnsigned);
	syt.addFlag("be", "bothExtend", MSyntax::kDouble);
	MArgParser argp = MArgParser(syt, arglist, &stat);
	if (argp.isFlagSet("curveSpans"))	stat = argp.getFlagArgument("curveSpans", 0, crvptnum);
	if (argp.isFlagSet("bothExtend"))	stat = argp.getFlagArgument("bothExtend", 0, bextend);
	return redoIt();
}

MStatus YourselfCommand::redoIt()
{
	MStatus suc;
	MItSelectionList itselections(selections);
	for (; itselections.isDone() != true; itselections.next()) {
		MDagPath tdag;
		itselections.getDagPath(tdag);
		//首先获取延展轴向和相应两个平面空间上的点
		MDoubleArray planepts12, planepts1, planepts2;//两个平面数据
		unsigned int spreadAxies;//延展轴向
		MDoubleArray mmvs;//延展轴向的最大距离
		suc = Get2DPoints(tdag, spreadAxies, planepts12, mmvs);
		if (!suc)MGlobal::displayInfo("--------------- Get2DPoints Error");
		bool toscale = false;
		MFnTransform fnTrsform(tdag);
		if ((mmvs[1] - mmvs[0]) < 10) {//如果数据太小,由于计算精度限制效果误差会较大,那么就以世界原点为中心放大100倍,计算完成再缩小回来
			double scales[3] = { 100 ,100 ,100 };
			fnTrsform.scaleBy(scales);
			toscale = true;
			suc = Get2DPoints(tdag, spreadAxies, planepts12, mmvs);
		}
		suc = HalfSplit(planepts12, planepts1, planepts2);
		Z5Matrix AM1 = Z5Matrix();
		Z5Matrix AM2 = Z5Matrix();
		Z5Vector YV1 = Z5Vector();
		Z5Vector YV2 = Z5Vector();
		suc = GetAugmentedMatrix(planepts1, AM1, YV1);
		suc = GetAugmentedMatrix(planepts2, AM2, YV2);
		Z5Vector Coe1 = AM1.InverseMatrix() * YV1;
		Z5Vector Coe2 = AM2.InverseMatrix() * YV2;
		MPointArray eps;
		//MGlobal::displayInfo(MString("mmin ") + mmvs[0]);
		//MGlobal::displayInfo(MString("mmax ") + mmvs[1]);
		suc = GenerateCrvEPs(Coe1, Coe2, spreadAxies, mmvs, crvptnum, eps);
		MFnNurbsCurve crvFn;
		MObject fcur = crvFn.createWithEditPoints(eps, 3, MFnNurbsCurve::kOpen, false, true, true);
		MFnDagNode tnode(crvFn.parent(0));
		tnode.setName(tdag.partialPathName() + "_fitcrv");
		if (toscale) {//如果数据太小,由于计算精度限制效果误差会较大,那么就以世界原点为中心放大100倍,计算完成再缩小回来
			MFnTransform tFnTrsform(crvFn.parent(0));
			double scales[3] = { 0.01, 0.01, 0.01 };
			tFnTrsform.setScalePivot(fnTrsform.scalePivot(MSpace::kTransform), MSpace::kTransform, true);
			tFnTrsform.scaleBy(scales);
			fnTrsform.scaleBy(scales);
		}
	}


	return MS::kSuccess;
}

MStatus YourselfCommand::undoIt()
{
	MGlobal::displayInfo("EmptyP command undone!\n");

	return MS::kSuccess;
}

bool YourselfCommand::isUndoable() const
{
	return true;
}

void* YourselfCommand::creator()
{
	return new YourselfCommand();
}

MStatus YourselfCommand::Get2DPoints(const MDagPath& inpath, unsigned int& saxies, MDoubleArray& oarray, MDoubleArray& mmv)
{
	MStatus succ;
	MFnMesh tfnmesh(inpath);
	MPointArray allpoints;
	succ = tfnmesh.getPoints(allpoints, MSpace::kWorld);
	unsigned int ptnum = tfnmesh.numVertices();
	if (ptnum < 3)	{
		MGlobal::displayError("Too few points");
		return MStatus::kFailure;
	}
	MDoubleArray xdatas ;
	MDoubleArray ydatas ;
	MDoubleArray zdatas ;
	MDoubleArray xzdatas = MDoubleArray(2*ptnum,0);//x轴
	MDoubleArray xydatas = MDoubleArray(2*ptnum,0);//x轴
	MDoubleArray yxdatas = MDoubleArray(2*ptnum,0);//y轴
	MDoubleArray yzdatas = MDoubleArray(2*ptnum,0);//y轴
	MDoubleArray zxdatas = MDoubleArray(2*ptnum,0);//z轴
	MDoubleArray zydatas = MDoubleArray(2*ptnum,0);//z轴
	for (int i = 0; i < ptnum; ++i) {
		succ = SortedAppend(allpoints[i].x, xdatas);
		succ = SortedAppend(allpoints[i].y, ydatas);
		succ = SortedAppend(allpoints[i].z, zdatas);
		xydatas[i * 2] = allpoints[i].x;  xydatas[i * 2 + 1] = allpoints[i].y;
		xzdatas[i * 2] = allpoints[i].x;  xzdatas[i * 2 + 1] = allpoints[i].z;
		yxdatas[i * 2] = allpoints[i].y;  yxdatas[i * 2 + 1] = allpoints[i].x;
		yzdatas[i * 2] = allpoints[i].y;  yzdatas[i * 2 + 1] = allpoints[i].z;
		zxdatas[i * 2] = allpoints[i].z;  zxdatas[i * 2 + 1] = allpoints[i].x;
		zydatas[i * 2] = allpoints[i].z;  zydatas[i * 2 + 1] = allpoints[i].y;
	}
	double xmm = xdatas[ptnum - 1] - xdatas[0];
	double ymm = ydatas[ptnum - 1] - ydatas[0];
	double zmm = zdatas[ptnum - 1] - zdatas[0];
	mmv.setLength(2);
	if (xmm > ymm) {
		if (xmm > zmm) {
			saxies = 0;
			oarray = xydatas + xzdatas;
			mmv[0] = xdatas[0];
			mmv[1] = xdatas[ptnum-1];
			return MStatus::kSuccess;
		}
		else {
			saxies = 2;
			oarray = zxdatas + zydatas;
			mmv[0] = zdatas[0];
			mmv[1] = zdatas[ptnum - 1];
			return MStatus::kSuccess;
		}
	}
	else if(ymm>zmm){
		saxies = 1;
		oarray = yxdatas + yzdatas;
		mmv[0] = ydatas[0];
		mmv[1] = ydatas[ptnum - 1];
		return MStatus::kSuccess;
	}
	else {
		saxies = 2;
		oarray = zxdatas + zydatas;
		mmv[0] = zdatas[0];
		mmv[1] = zdatas[ptnum - 1];
		return MStatus::kSuccess;
	}
	return succ;
}

MStatus YourselfCommand::SortedAppend(double invalue, MDoubleArray& darray)
{
	if (darray.length() < 1)darray.append(invalue);
	else {
		for (int j = 0; j < darray.length(); ++j) {
			if (invalue < darray[j]) {
				darray.insert(invalue, j);
				break;
			}
			if (j == darray.length() - 1){
				darray.append(invalue);
				break;
			}
		}
	}
	return MStatus::kSuccess;
}

MStatus YourselfCommand::GetAugmentedMatrix(MDoubleArray inpts, Z5Matrix& oMat, Z5Vector& oVec)
{
	double S_xi0 = 0;
	double S_xi1 = 0;
	double S_xi2 = 0;
	double S_xi3 = 0;
	double S_xi4 = 0;
	double S_xi5 = 0;
	double S_xi6 = 0;
	double S_xi7 = 0;
	double S_xi8 = 0;
	MDoubleArray S_yi = MDoubleArray(5, 0);
	unsigned int halflen = inpts.length() / 2;
	for (int i = 0; i < halflen; ++i) {
		S_xi0 += 1;
		S_xi1 += inpts[i * 2];
		S_xi2 += pow(inpts[i * 2], 2);
		S_xi3 += pow(inpts[i * 2], 3);
		S_xi4 += pow(inpts[i * 2], 4);
		S_xi5 += pow(inpts[i * 2], 5);
		S_xi6 += pow(inpts[i * 2], 6);
		S_xi7 += pow(inpts[i * 2], 7);
		S_xi8 += pow(inpts[i * 2], 8);
		S_yi[0] += inpts[i * 2 + 1];
		S_yi[1] += (inpts[i * 2] * inpts[i * 2 + 1]);
		S_yi[2] += (pow(inpts[i * 2], 2) * inpts[i * 2 + 1]);
		S_yi[3] += (pow(inpts[i * 2], 3) * inpts[i * 2 + 1]);
		S_yi[4] += (pow(inpts[i * 2], 4) * inpts[i * 2 + 1]);
	}
		oMat.SetElement(0, 0, S_xi0);
		oMat.SetElement(0, 1, S_xi1);
		oMat.SetElement(0, 2, S_xi2);
		oMat.SetElement(0, 3, S_xi3);
		oMat.SetElement(0, 4, S_xi4);
		oMat.SetElement(1, 0, S_xi1);
		oMat.SetElement(1, 1, S_xi2);
		oMat.SetElement(1, 2, S_xi3);
		oMat.SetElement(1, 3, S_xi4);
		oMat.SetElement(1, 4, S_xi5);
		oMat.SetElement(2, 0, S_xi2);
		oMat.SetElement(2, 1, S_xi3);
		oMat.SetElement(2, 2, S_xi4);
		oMat.SetElement(2, 3, S_xi5);
		oMat.SetElement(2, 4, S_xi6);
		oMat.SetElement(3, 0, S_xi3);
		oMat.SetElement(3, 1, S_xi4);
		oMat.SetElement(3, 2, S_xi5);
		oMat.SetElement(3, 3, S_xi6);
		oMat.SetElement(3, 4, S_xi7);
		oMat.SetElement(4, 0, S_xi4);
		oMat.SetElement(4, 1, S_xi5);
		oMat.SetElement(4, 2, S_xi6);
		oMat.SetElement(4, 3, S_xi7);
		oMat.SetElement(4, 4, S_xi8);

		oVec.SetElements(S_yi);
	return MStatus::kSuccess;
}

MStatus YourselfCommand::HalfSplit(MDoubleArray inpts12, MDoubleArray& inpts1, MDoubleArray& inpts2)
{
	unsigned int halflen = inpts12.length() / 2;
	inpts1.setLength(halflen);
	inpts2.setLength(halflen);
	for (int i = 0; i < halflen; ++i) {
		inpts1[i] = inpts12[i];
		inpts2[i] = inpts12[i + halflen];
	}
	return MStatus::kSuccess;
}

MStatus YourselfCommand::GenerateCrvEPs(Z5Vector inc1, Z5Vector inc2, unsigned int spreadAxies, const MDoubleArray& mmv, unsigned int epnum, MPointArray& oPoints)
{
	oPoints.setLength(epnum + 1);
	//延展轴向前后多延展5%
	double minv = mmv[0] - (mmv[1] - mmv[0]) * bextend;
	double maxv = mmv[1] + (mmv[1] - mmv[0]) * bextend;
	double seglength = (maxv - minv) / epnum;
	switch (spreadAxies)
	{
	case 0:
		for (unsigned int i = 0; i < epnum + 1; ++i) {
			double fi = minv + i * seglength;
			double tv1 = inc1.GetElement(0) + inc1.GetElement(1) * fi + inc1.GetElement(2) * pow(fi, 2) + inc1.GetElement(3) * pow(fi, 3) + inc1.GetElement(4) * pow(fi, 4);
			double tv2 = inc2.GetElement(0) + inc2.GetElement(1) * fi + inc2.GetElement(2) * pow(fi, 2) + inc2.GetElement(3) * pow(fi, 3) + inc2.GetElement(4) * pow(fi, 4);
			oPoints[i] = MPoint(fi, tv1, tv2);
		}
		break;
	case 1:
		for (unsigned int i = 0; i < epnum + 1; ++i) {
			double fi = minv + i * seglength;
			double tv1 = inc1.GetElement(0) + inc1.GetElement(1) * fi + inc1.GetElement(2) * pow(fi, 2) + inc1.GetElement(3) * pow(fi, 3) + inc1.GetElement(4) * pow(fi, 4);
			double tv2 = inc2.GetElement(0) + inc2.GetElement(1) * fi + inc2.GetElement(2) * pow(fi, 2) + inc2.GetElement(3) * pow(fi, 3) + inc2.GetElement(4) * pow(fi, 4);
			oPoints[i] = MPoint(tv1, fi, tv2);
		}
		break;
	case 2:
		for (unsigned int i = 0; i < epnum + 1; ++i) {
			double fi = minv + i * seglength;
			double tv1 = inc1.GetElement(0) + inc1.GetElement(1) * fi + inc1.GetElement(2) * pow(fi, 2) + inc1.GetElement(3) * pow(fi, 3) + inc1.GetElement(4) * pow(fi, 4);
			double tv2 = inc2.GetElement(0) + inc2.GetElement(1) * fi + inc2.GetElement(2) * pow(fi, 2) + inc2.GetElement(3) * pow(fi, 3) + inc2.GetElement(4) * pow(fi, 4);
			oPoints[i] = MPoint(tv1, tv2, fi);
		}
		break;
	default:
		break;
	}
	return MStatus::kSuccess;
}
