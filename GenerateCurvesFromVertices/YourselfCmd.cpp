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
	if (argp.isFlagSet("rowcurveSpans"))	stat = argp.getFlagArgument("rowcurveSpans", 0, rowcurveSpans);
	if (argp.isFlagSet("bothExtend"))	stat = argp.getFlagArgument("bothExtend", 0, bextend);
	return redoIt();
}

MStatus YourselfCommand::redoIt()
{
	MStatus suc;
	MItSelectionList itselections(selections);
	MStringArray results;
	for (; itselections.isDone() != true; itselections.next()) {
		MDagPath tdag;
		itselections.getDagPath(tdag);
		MFnMesh tfnmesh(tdag);
		
		MPointArray allpoints;
		suc = tfnmesh.getPoints(allpoints, MSpace::kWorld);
		MPointArray maincrvEPs = GenerateMatchedCurveEPs(allpoints,false, crvptnum);
		MPointArray mainUniEPs;
		MDagPath crvpath = GenerateEPCurve(maincrvEPs, tdag.partialPathName() + "_fitcrv_M", mainUniEPs);
		
		MFnNurbsCurve crvFn(crvpath);
		MIntArray* splitedIndexes = SlitedPointsByCrv(allpoints, crvFn);
		//返回每一段中间位置的切线,正交化后的法线,和此切线法线组成的局部坐标系变换到世界坐标系的旋转四元数
		MVectorArray midpos,midtangents, midnormals;
		for (int i = 0; i < crvptnum; ++i) {
			//MGlobal::displayInfo(MString("splitedIndexes length ") + splitedIndexes[i].length());
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
		MDagPathArray parCrvs;
		MPointArray* uniEPss = new MPointArray[crvptnum];
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
			MPointArray partCrvEPs = GenerateMatchedCurveEPs(partpoints, true, rowcurveSpans);
			for (MPoint& ep : partCrvEPs) {
				ep *= qf.inverse();
			}
			MPointArray tunieps;
			parCrvs.append(GenerateEPCurve(partCrvEPs, tdag.partialPathName() + MString("_partfitcrv") + i, tunieps));
			uniEPss[i] = tunieps;
		}
		MDagPath simmeshpath = ZGenMesh(GenMeshParam(uniEPss, crvptnum), tdag.partialPathName() + MString("_fitmesh"));
		results.append(simmeshpath.partialPathName());
		//删除辅助曲线
		MObject maincrvobj = crvpath.node();
		MGlobal::deleteNode(maincrvobj);
		for (MDagPath tpath : parCrvs) {
			MObject tobj= tpath.node();
			MGlobal::deleteNode(tobj);
		}
	}
	setResult(results);
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

MStatus YourselfCommand::FitMesh(const MPointArray& allpoints,bool makeXsp, MPointArray& outpoints, MPoint& meshCenter, double& scalefactor,unsigned int& spAxis,MDoubleArray& mmvs)
{
	unsigned int ptnum = allpoints.length();
	if (ptnum < 3) {
		MGlobal::displayError("Too few points");
		return MStatus::kFailure;
	}
	MDoubleArray xdatas;
	MDoubleArray ydatas;
	MDoubleArray zdatas;
	for (MPoint tpt : allpoints) {
		SortedAppend(tpt.x, xdatas);
		SortedAppend(tpt.y, ydatas);
		SortedAppend(tpt.z, zdatas);
	}
	meshCenter.x = 0.5 * (xdatas[ptnum - 1] + xdatas[0]);
	meshCenter.y = 0.5 * (ydatas[ptnum - 1] + ydatas[0]);
	meshCenter.z = 0.5 * (zdatas[ptnum - 1] + zdatas[0]);
	if (makeXsp) {
		double xmm = xdatas[ptnum - 1] - xdatas[0];
		spAxis = 0;
		scalefactor = 1;
		while (xmm < 10)
		{
			xmm *= 100; scalefactor *= 100;
		}
		while (xmm > 1000)
		{
			xmm *= 0.1; scalefactor *= 0.1;
		}
		mmvs.setLength(2);
		mmvs[0] = scalefactor * (xdatas[0] - meshCenter.x);
		mmvs[1] = scalefactor * (xdatas[ptnum - 1] - meshCenter.x);

	}
	else {
		double xmm = xdatas[ptnum - 1] - xdatas[0];
		double ymm = ydatas[ptnum - 1] - ydatas[0];
		double zmm = zdatas[ptnum - 1] - zdatas[0];
		if (xmm > ymm) {
			if (xmm > zmm) {
				spAxis = 0;
				scalefactor = 1;
				while (xmm < 10)
				{
					xmm *= 100; scalefactor *= 100;
				}
				while (xmm > 1000)
				{
					xmm *= 0.1; scalefactor *= 0.1;
				}
				mmvs.setLength(2);
				mmvs[0] = scalefactor * (xdatas[0] - meshCenter.x);
				mmvs[1] = scalefactor * (xdatas[ptnum - 1] - meshCenter.x);
			}
			else {
				spAxis = 2;
				scalefactor = 1;
				while (zmm < 10)
				{
					zmm *= 100; scalefactor *= 100;
				}
				while (zmm > 1000)
				{
					zmm *= 0.1; scalefactor *= 0.1;
				}
				mmvs.setLength(2);
				mmvs[0] = scalefactor * (zdatas[0] - meshCenter.z);
				mmvs[1] = scalefactor * (zdatas[ptnum - 1] - meshCenter.z);
			}
		}
		else if (ymm > zmm) {
			spAxis = 1;
			scalefactor = 1;
			while (ymm < 10)
			{
				ymm *= 100; scalefactor *= 100;
			}
			while (ymm > 1000)
			{
				ymm *= 0.1; scalefactor *= 0.1;
			}
			mmvs.setLength(2);
			mmvs[0] = scalefactor * (ydatas[0] - meshCenter.y);
			mmvs[1] = scalefactor * (ydatas[ptnum - 1] - meshCenter.y);
		}
		else {
			spAxis = 2;
			scalefactor = 1;
			while (zmm < 10)
			{
				zmm *= 100; scalefactor *= 100;
			}
			while (zmm > 1000)
			{
				zmm *= 0.1; scalefactor *= 0.1;
			}
			mmvs.setLength(2);
			mmvs[0] = scalefactor * (zdatas[0] - meshCenter.z);
			mmvs[1] = scalefactor * (zdatas[ptnum - 1] - meshCenter.z);
		}
	}
	
	for (MPoint pt : allpoints) {
		outpoints.append((pt - meshCenter) * scalefactor);//先移动再缩放
	}
	return MStatus::kSuccess;
}

MStatus YourselfCommand::Get2DPoints(const MPointArray& allpoints, const unsigned int& saxies, MDoubleArray& oarray1, MDoubleArray& oarray2)
{
	unsigned int ptnum = allpoints.length();
	oarray1.setLength(2*ptnum); oarray2.setLength(2*ptnum);
	switch (saxies)
	{
	case 0:
		for (int i = 0; i < ptnum; ++i) {
			oarray1[i * 2] = allpoints[i].x;  oarray1[i * 2 + 1] = allpoints[i].y;
			oarray2[i * 2] = allpoints[i].x;  oarray2[i * 2 + 1] = allpoints[i].z;
		}
		break;
	case 1:
		for (int i = 0; i < ptnum; ++i) {
			oarray1[i * 2] = allpoints[i].y;  oarray1[i * 2 + 1] = allpoints[i].x;
			oarray2[i * 2] = allpoints[i].y;  oarray2[i * 2 + 1] = allpoints[i].z;
		}
		break;
	case 2:
		for (int i = 0; i < ptnum; ++i) {
			oarray1[i * 2] = allpoints[i].z;  oarray1[i * 2 + 1] = allpoints[i].x;
			oarray2[i * 2] = allpoints[i].z;  oarray2[i * 2 + 1] = allpoints[i].y;
		}
		break;
	default:
		break;
	}
	return MStatus::kSuccess;
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

MPointArray YourselfCommand::GenerateCrvEPs(Z5Vector inc1, Z5Vector inc2, unsigned int spreadAxies, const MDoubleArray& mmv, unsigned int epnum)
{
	MPointArray oPoints;
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
	return oPoints;
}

MIntArray* YourselfCommand::SlitedPointsByCrv(const MPointArray& allpoints, const MFnNurbsCurve& incrv)
{
	MIntArray* tsplitedIndexes = new MIntArray[crvptnum];//这里可以从incrv里获取其段数
	for (int i = 0; i < allpoints.length(); ++i) {
		double param;
		MPoint cpt = incrv.closestPoint(allpoints[i], true, &param, kMFnNurbsEpsilon, MSpace::kObject, NULL);
		//tfnmesh.setVertexColor(colorsArray[param], i);
		//MGlobal::displayInfo(MString("param ") + param);
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

MPointArray YourselfCommand::GenerateMatchedCurveEPs(const MPointArray& allpoints,  bool makeXsp, unsigned int epcount)
{
	if (allpoints.length() < 3) {
		MGlobal::displayError("Too few points");
		return MPointArray();
	}
	//首先获取延展轴向和相应两个平面空间上的点
	MDoubleArray planepts1, planepts2;//两个平面数据
	MDoubleArray mmvs;//延展轴向的最大距离
	MPoint meshCenter;
	MPointArray fittedpoints;
	double scalefactor;
	unsigned int spAxises;
	FitMesh(allpoints, makeXsp, fittedpoints, meshCenter, scalefactor, spAxises, mmvs);
	Get2DPoints(fittedpoints, spAxises, planepts1, planepts2);
	Z5Matrix AM1 = Z5Matrix();
	Z5Matrix AM2 = Z5Matrix();
	Z5Vector YV1 = Z5Vector();
	Z5Vector YV2 = Z5Vector();
	GetAugmentedMatrix(planepts1, AM1, YV1);
	GetAugmentedMatrix(planepts2, AM2, YV2);
	Z5Vector Coe1 = AM1.InverseMatrix() * YV1;
	Z5Vector Coe2 = AM2.InverseMatrix() * YV2;
	MPointArray eps = GenerateCrvEPs(Coe1, Coe2, spAxises, mmvs, epcount);
	double scalefactor_1 = 1.0 / scalefactor;
	for (MPoint& tep : eps) {//先缩放后移动
		if (makeXsp) {//沿y轴向向中心点收缩50%
			tep.y += (meshCenter.y - tep.y)*0.5;
		}
		tep *= scalefactor_1;
		tep += meshCenter;
	}
	return eps;
}

MDagPath YourselfCommand::GenerateEPCurve(const MPointArray& ineps, MString inname,MPointArray& unieps)
{
	unsigned int spans = ineps.length()-1;
	MFnNurbsCurve partcrvFn;
	MObject tcrv = partcrvFn.createWithEditPoints(ineps, 3, MFnNurbsCurve::kOpen, false, true, true);
	double epartlen = partcrvFn.length()/spans;//线段总长度/段数(每个分段的平均长度)
	for (int i = 0; i < ineps.length(); ++i) {
		double tpartparam = partcrvFn.findParamFromLength(i * epartlen);
		MPoint newEp;
		partcrvFn.getPointAtParam(tpartparam, newEp, MSpace::kWorld);
		unieps.append(newEp);
	}
	MGlobal::deleteNode(tcrv);
	partcrvFn.createWithEditPoints(unieps, 3, MFnNurbsCurve::kOpen, false, true, true);//重建均匀分布的EP曲线
	MFnDagNode tpartnode(partcrvFn.parent(0));
	tpartnode.setName(inname);
	MDagPath tdagp;
	tpartnode.getPath(tdagp);
	return tdagp;
}

ZGenMeshParam YourselfCommand::GenMeshParam(MPointArray* inpointArrays,unsigned int ptarrayNum)
{
	ZGenMeshParam fpara = ZGenMeshParam();
	unsigned int vecols = inpointArrays[0].length();
	fpara.numVertices = ptarrayNum * vecols;
	unsigned int pgrows = ptarrayNum - 1;
	unsigned int pgcols = vecols - 1;
	fpara.numPolygons = pgrows * pgcols;
	fpara.polygonCounts = MIntArray(fpara.numPolygons, 4);
	fpara.vertexArray.setLength(fpara.numVertices);
	for (int i = 0; i < ptarrayNum; ++i) {
		for (int j = 0; j < vecols; ++j) {
			fpara.vertexArray[i * vecols + j] = inpointArrays[i][j];
		}
	}
	fpara.polygonConnects.setLength(4 * fpara.numPolygons);
	for (int i = 0; i < pgrows; ++i) {
		for (int j = 0; j < pgcols; ++j) {
			unsigned int faceidx = 4 * (pgcols * i + j);
			fpara.polygonConnects[faceidx] = vecols * i + j;
			fpara.polygonConnects[faceidx + 1] = vecols * i + j + 1;
			fpara.polygonConnects[faceidx + 2] = vecols * (i + 1) + j + 1;
			fpara.polygonConnects[faceidx + 3] = vecols * (i + 1) + j;
		}
	}
	return fpara;
}

MDagPath YourselfCommand::ZGenMesh(const ZGenMeshParam& inpa,MString inname)
{
	MFnMesh tfnmesh;
	tfnmesh.create(inpa.numVertices,
										inpa.numPolygons, 
										inpa.vertexArray, 
										inpa.polygonCounts, 
										inpa.polygonConnects, 
										inpa.parentOrOwner, 
										inpa.ReturnStatus);
	MFnDagNode tpartnode(tfnmesh.parent(0));
	tpartnode.setName(inname);
	MDagPath tdagp;
	tpartnode.getPath(tdagp);
	return tdagp;
	//MGlobal::displayInfo(meshpath.fullPathName());
	//MMaterial defaultmat = MMaterial::defaultMaterial();
	//defaultmat.setMaterial(meshpath,true);
}

