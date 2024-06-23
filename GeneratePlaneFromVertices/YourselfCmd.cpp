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
	syt.addFlag("sp", "curveSpans", MSyntax::kUnsigned);
	syt.addFlag("rsp", "rowcurveSpans", MSyntax::kUnsigned);
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
		MFnMesh tfnmesh(tdag);
		MPointArray allpoints;
		suc = tfnmesh.getPoints(allpoints, MSpace::kWorld);
		//首先获取延展轴向和相应两个平面空间上的点
		MDoubleArray planepts12, planepts1, planepts2;//两个平面数据
		unsigned int spreadAxies;//延展轴向
		MDoubleArray mmvs;//延展轴向的最大距离
		suc = Get2DPoints(allpoints, spreadAxies, planepts12, mmvs);
		if (!suc)MGlobal::displayInfo("--------------- Get2DPoints Error");
		bool toscale = false;
		MFnTransform fnTrsform(tdag);
		if ((mmvs[1] - mmvs[0]) < 10) {//如果数据太小,由于计算精度限制效果误差会较大,那么就以世界原点为中心放大100倍,计算完成再缩小回来
			double scales[3] = { 100 ,100 ,100 };
			fnTrsform.scaleBy(scales);
			toscale = true;
			tfnmesh.getPoints(allpoints, MSpace::kWorld);
			suc = Get2DPoints(allpoints, spreadAxies, planepts12, mmvs);
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
		//MColorArray colorsArray;
		//for (int i = 0; i < crvptnum+1; ++i) {
		//	float	 tcs[3];
		//	MRandom::Rand_3f(tcs, i, 9999);
		//	colorsArray.append(MColor(tcs));
		//}
		MIntArray* splitedIndexes = SlitedPointsByCrv(allpoints, crvFn);
		//返回每一段中间位置的切线,正交化后的法线,和此切线法线组成的局部坐标系变换到世界坐标系的旋转四元数
		MVectorArray midpos,midtangents, midnormals;
		for (int i = 0; i < crvptnum; ++i) {
			MPoint tmpt;
			crvFn.getPointAtParam(i + 0.5, tmpt);
			midpos.append(tmpt);
			MVector tnormal = crvFn.normal(i + 0.5).normal();
			MVector ttangent = crvFn.tangent(i + 0.5).normal();
			double tmaxdistance = 0.0;
			int maxdisInd = 0;
			MPoint closetP;
			for (int ind : splitedIndexes[i]) {
				MPoint tclosetP = crvFn.closestPoint(allpoints[ind]);
				double tdis = (allpoints[ind] - tclosetP).length();
				if (tdis > tmaxdistance) {
					tmaxdistance = tdis;
					maxdisInd = ind;
					closetP = tclosetP;
				}
			}
			MVector desv = (allpoints[maxdisInd] - closetP).normal();
			desv = (ttangent ^ desv) ^ ttangent;//施密特正交化
			tnormal = tnormal.rotateBy(tnormal.rotateTo(desv));
			midnormals.append(tnormal);
			midtangents.append(ttangent);
		}
		StraightenNormals(midnormals, midtangents);
		MQuaternion* toGquas = new MQuaternion[crvptnum];
		for (int i = 0; i < crvptnum; ++i) {
			//求局部坐标系到全局坐标系的旋转四元数
			MQuaternion tq1 = midnormals[i].rotateTo(MVector(1, 0, 0));
			MVector tttangent = midtangents[i] * tq1;
			MQuaternion tq2 = tttangent.rotateTo(MVector(0, 1, 0));
			MQuaternion qf = tq1 * tq2;

			unsigned int partlen = splitedIndexes[i].length();
			MPointArray partpoints(partlen);
			for (int j = 0; j < partlen; ++j) {
				partpoints[j] = allpoints[splitedIndexes[i][j]] * qf;
			}
			MPointArray EPs=GeneratePartCrvEPs(partpoints);
			for (MPoint& ep : EPs) {
				ep *= qf.inverse();
			}
			MFnNurbsCurve partcrvFn;
			MObject fcur = partcrvFn.createWithEditPoints(EPs, 3, MFnNurbsCurve::kOpen, false, true, true);
			MFnDagNode tpartnode(partcrvFn.parent(0));
			tpartnode.setName(tdag.partialPathName() + MString("_partfitcrv") + i);
		}
		for (int i = 0; i < crvptnum; ++i) {
			MGlobal::displayInfo(MString("[[") + midpos[i].x + "," + midpos[i].y + "," + midpos[i].z + "," + "],[" + midtangents[i].x + "," + midtangents[i].y + "," + midtangents[i].z + "],[" + midnormals[i].x + "," + midnormals[i].y + "," + midnormals[i].z + "]]");
		}
		//for (int i = 0; i < crvptnum; ++i) {
		//	for (int j = 0; j < splitedIndexes[i].length(); ++j) {
		//		tfnmesh.setVertexColor(colorsArray[i], splitedIndexes[i][j]);
		//	}
		//}
		//for (double item : allparams)MGlobal::displayInfo(MString("") + item);
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

MStatus YourselfCommand::Get2DPoints(const MPointArray& allpoints, unsigned int& saxies, MDoubleArray& oarray, MDoubleArray& mmv)
{
	MStatus succ;
	unsigned int ptnum = allpoints.length();
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

MIntArray* YourselfCommand::SlitedPointsByCrv(const MPointArray& allpoints, const MFnNurbsCurve& incrv)
{
	MIntArray* tsplitedIndexes = new MIntArray[crvptnum];//这里可以从incrv里获取其段数
	for (int i = 0; i < allpoints.length(); ++i) {
		double param;
		MPoint cpt = incrv.closestPoint(allpoints[i], true, &param, kMFnNurbsEpsilon, MSpace::kObject, NULL);
		//tfnmesh.setVertexColor(colorsArray[param], i);
		tsplitedIndexes[(int)param].append(i);
	}
	return tsplitedIndexes;
}

MStatus YourselfCommand::StraightenNormals(MVectorArray& innormals, const MVectorArray& intangents)
{
	unsigned int innlen = innormals.length();
	if (innlen > 2) {
		for (int i = 1; i < innlen; ++i) {
			if (innormals[i]*innormals[i - 1] < 0) {
				innormals[i] = innormals[i].rotateBy(MQuaternion(3.1415926, intangents[i]));
			}
		}
		return MStatus::kSuccess;
	}
	return MStatus::kFailure;
}

MPointArray YourselfCommand::GeneratePartCrvEPs(const MPointArray& inpoints)
{
	if (inpoints.length() < 2) { MGlobal::displayError("分段点数太少!"); return MPointArray(); }
	int ptnum = inpoints.length();
	MDoubleArray xydatas(2 * ptnum, 0.0), xzdatas(2 * ptnum, 0.0), xdatas;
	for (int i = 0; i < ptnum; ++i) {
		SortedAppend(inpoints[i].x, xdatas);
		xydatas[i * 2] = inpoints[i].x; xydatas[i * 2 + 1] = inpoints[i].y;
		xzdatas[i * 2] = inpoints[i].x; xzdatas[i * 2 + 1] = inpoints[i].z;
	}
	MDoubleArray mmvs;
	mmvs.append(xdatas[0]); mmvs.append(xdatas[ptnum - 1]);
	bool bescale = false;
	if (xdatas[ptnum - 1] - xdatas[0] < 10) {
		for (int i = 0; i < 2 * ptnum; ++i) {
			xydatas[i] *= 100; xzdatas[i] *= 100;
		}
		mmvs[0] *= 100; mmvs[1] *= 100;
		bescale = true;
		//MGlobal::displayInfo("bescale true!");
	}
	Z5Matrix AM1 = Z5Matrix();
	Z5Matrix AM2 = Z5Matrix();
	Z5Vector YV1 = Z5Vector();
	Z5Vector YV2 = Z5Vector();
	GetAugmentedMatrix(xydatas, AM1, YV1);
	GetAugmentedMatrix(xzdatas, AM2, YV2);
	Z5Vector Coe1 = AM1.InverseMatrix() * YV1;
	Z5Vector Coe2 = AM2.InverseMatrix() * YV2;
	MPointArray teps;
	GenerateCrvEPs(Coe1, Coe2, 0, mmvs, rowcurveSpans,teps);
	if (bescale) {
		for (MPoint& ep : teps) {
			ep *= 0.01;
		}
	}
	return teps;
}
