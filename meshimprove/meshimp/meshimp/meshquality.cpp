#include"meshquality.h"
#include<cmath>
#include<iostream>
using namespace std;
#define PI 3.14159265

vector<double> computesizelenQuality(MyMesh& mesh)
{
	vector<double> qualities;
	for (auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
	{
		MyMesh::Point p0;
		MyMesh::Point p1;
		MyMesh::Point p2;
		auto fv_it = mesh.fv_iter(*f_it);
		p0 = mesh.point(*fv_it);
		++fv_it;
		p1 = mesh.point(*fv_it);
		++fv_it;
		p2 = mesh.point(*fv_it);
		MyMesh::Normal v0 = p1 - p0;
		MyMesh::Normal v1 = p2 - p0;
		MyMesh::Normal v2 = p1 - p2;
		MyMesh::Scalar area = (v0%v1).length();
		MyMesh::Scalar quality = (v0.length()*v0.length() + v1.length()*v1.length() + v2.length()*v2.length())/( area );
		qualities.push_back(quality);
	}
	return qualities;
}
vector<double> computeMeshIdealElementQuality(MyMesh& mesh)
{
	vector<double> qualities;
	for (auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
	{
		MyMesh::Point p0;
		MyMesh::Point p1;
		MyMesh::Point p2;
		auto fv_it = mesh.fv_iter(*f_it);
		p0 = mesh.point(*fv_it);
		++fv_it;
		p1 = mesh.point(*fv_it);
		++fv_it;
		p2 = mesh.point(*fv_it);
		MyMesh::Normal v0 = p1-p0;
		MyMesh::Normal v1 = p2-p0;
		MyMesh::Normal v2 = p1-p2;
		MyMesh::Scalar area = (v0%v1).length();
		MyMesh::Scalar quality = 2 * sqrt(3) * area / (v0.length()*v0.length() + v1.length()*v1.length() + v2.length()*v2.length());
		qualities.push_back(quality);
	}
	return qualities;
}

vector<double> computeMeshIdealElementQuality(MyTetMesh& mesh)
{
	vector<double> qualities;
	for (auto c_it = mesh.cells_begin(); c_it != mesh.cells_end(); c_it++)
	{
		Vec3d p0;
		Vec3d p1;
		Vec3d p2;
		Vec3d p3;
		auto cv_it = mesh.cv_iter(*c_it);
		p0 = mesh.vertex(*cv_it);
		++cv_it;
		p1 = mesh.vertex(*cv_it);
		++cv_it;
		p2 = mesh.vertex(*cv_it);
		++cv_it;
		p3 = mesh.vertex(*cv_it);
		Vec3d v0 = p1 - p0;
		Vec3d v1 = p2 - p0;
		Vec3d v2 = p1 - p2;
		Vec3d v3 = p3 - p0;
		Vec3d v4 = p3 - p1;
		Vec3d v5 = p3 - p2;
		double volume = abs(((v0%v1) | v3));
		double quality = 12 * sqrt(3) * volume / pow((v0.length()*v0.length() + v1.length()*v1.length() + v2.length()*v2.length() + v3.length()*v3.length() + v4.length()*v4.length() + v5.length()*v5.length()),1.5);
		qualities.push_back(quality);
	}
	return qualities;
}


vector<double> computeMeshAngleQuality(MyMesh& mesh)
{
	vector<double> qualities;
	for (auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
	{
		MyMesh::Point p0;
		MyMesh::Point p1;
		MyMesh::Point p2;
		MyMesh::Normal v0;
		MyMesh::Normal v1;
		MyMesh::Scalar cos2a;
		MyMesh::Scalar angle;
		auto fv_it = mesh.fv_iter(*f_it);
		p0 = mesh.point(*fv_it);
		++fv_it;
		p1 = mesh.point(*fv_it);
		++fv_it;
		p2 = mesh.point(*fv_it);
		
		v0 = p1 - p0;
		v1 = p2 - p0;
		v0.normalize();
		v1.normalize();
		cos2a = v0 | v1;
		angle = acos(cos2a)*180. / PI;
		qualities.push_back(angle);

		v0 = p0 - p1;
		v1 = p2 - p1;
		v0.normalize();
		v1.normalize();
		cos2a = v0 | v1;
		angle = acos(cos2a)*180. / PI;
		qualities.push_back(angle);

		v0 = p1 - p2;
		v1 = p0 - p2;
		v0.normalize();
		v1.normalize();
		cos2a = v0 | v1;
		angle = acos(cos2a)*180. / PI;
		qualities.push_back(angle);
	}
	return qualities;
}

double computeMeshAveang(MyMesh& mesh)
{
	vector<double> qualities;
	for (auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
	{
		double minang = 1000.;
		MyMesh::Point p0;
		MyMesh::Point p1;
		MyMesh::Point p2;
		MyMesh::Normal v0;
		MyMesh::Normal v1;
		MyMesh::Scalar cos2a;
		MyMesh::Scalar angle;
		auto fv_it = mesh.fv_iter(*f_it);
		p0 = mesh.point(*fv_it);
		++fv_it;
		p1 = mesh.point(*fv_it);
		++fv_it;
		p2 = mesh.point(*fv_it);

		v0 = p1 - p0;
		v1 = p2 - p0;
		v0.normalize();
		v1.normalize();
		cos2a = v0 | v1;
		angle = acos(cos2a)*180. / PI;
		if (angle < minang)
			minang = angle;
		//qualities.push_back(angle);

		v0 = p0 - p1;
		v1 = p2 - p1;
		v0.normalize();
		v1.normalize();
		cos2a = v0 | v1;
		angle = acos(cos2a)*180. / PI;
		//qualities.push_back(angle);
		if (angle < minang)
			minang = angle;

		v0 = p1 - p2;
		v1 = p0 - p2;
		v0.normalize();
		v1.normalize();
		cos2a = v0 | v1;
		angle = acos(cos2a)*180. / PI;
		//qualities.push_back(angle);
		if (angle < minang)
			minang = angle;
		qualities.push_back(minang);
	}
	double aveang = 0;
	for (int i = 0; i < qualities.size(); i++)
	{
		aveang = aveang + qualities[i];
	}
	aveang = aveang / qualities.size();
	return aveang;
}


vector<double> computeMeshAngleQuality(PolyMesh& mesh)
{
	vector<double> qualities;
	for (auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
	{
		vector<OpenMesh::HalfedgeHandle> e_hs;
		auto heh_init = mesh.halfedge_handle(*f_it);
		e_hs.push_back(heh_init);
		auto heh = mesh.next_halfedge_handle(heh_init);
		while (heh != heh_init) {
			e_hs.push_back(heh);
			heh = mesh.next_halfedge_handle(heh);
		}
		PolyMesh::HalfedgeHandle e1;
		PolyMesh::HalfedgeHandle e2;
		for (int i = 0; i < e_hs.size(); i++)
		{
			e1 = e_hs[i];
			if (i == e_hs.size() - 1)
				e2 = e_hs[0];
			else
				e2 = e_hs[i + 1];

			auto p1 = mesh.from_vertex_handle(e1);
			auto p2 = mesh.to_vertex_handle(e1);
			auto p3 = mesh.to_vertex_handle(e2);

			OpenMesh::Vec3f v1 = (mesh.point(p1) - mesh.point(p2)).normalize();
			OpenMesh::Vec3f v2 = (mesh.point(p3) - mesh.point(p2)).normalize();

			MyMesh::Scalar cos2a = v1 | v2;
			MyMesh::Scalar angle = acos(cos2a)*180. / PI;

			qualities.push_back(angle);

		}
	}
	return qualities;
}



vector<double> computeMeshAngleQuality(MyTetMesh& mesh)
{
	vector<double> qualities;
	for (auto c_it = mesh.cells_begin(); c_it != mesh.cells_end(); c_it++)
	{
		vector<Vec3d> fnors;
		vector<OpenVolumeMesh::HalfFaceHandle> fhs = mesh.cell(*c_it).halffaces();
		//记录每个半面的法向量
		for (int i = 0; i < fhs.size(); i++)
		{
			vector<OpenVolumeMesh::VertexHandle> vhs = mesh.get_halfface_vertices(fhs[i]);
			Vec3d normal = ((mesh.vertex(vhs[1]) - mesh.vertex(vhs[0])) % (mesh.vertex(vhs[2]) - mesh.vertex(vhs[0]))).normalize();
			fnors.push_back(normal);
		}

		Vec3d normal0 = fnors[0];
		Vec3d normal1 = fnors[1];
		Vec3d normal2 = fnors[2];
		Vec3d normal3 = fnors[3];

		double cos2a;
		double angle;

		//第一个二面角
		cos2a = normal0 | normal1;
		angle = acos(-cos2a)*180. / PI;
		qualities.push_back(angle);
		//第二个二面角
		cos2a = normal0 | normal2;
		angle = acos(-cos2a)*180. / PI;
		qualities.push_back(angle);
		//第三个二面角
		cos2a = normal0 | normal3;
		angle = acos(-cos2a)*180. / PI;
		qualities.push_back(angle);
		//第四个二面角
		cos2a = normal1 | normal2;
		angle = acos(-cos2a)*180. / PI;
		qualities.push_back(angle);
		//第五个二面角
		cos2a = normal1 | normal2;
		angle = acos(-cos2a)*180. / PI;
		qualities.push_back(angle);
		//第六个二面角
		cos2a = normal2 | normal3;
		angle = acos(-cos2a)*180. / PI;
		qualities.push_back(angle);
	}
	return qualities;
}

double computeLocalAngleQuality(MyMesh::Point p, vector<MyMesh::Point> ring){
	double minangle = 360, maxangle = 0;
	for (int i = 0; i < ring.size(); i++){
		MyMesh::Normal v, v0, v1;
		v = ring[(i + 1) % ring.size()] - ring[i];
		v0 = p - ring[(i + 1) % ring.size()];
		v1 = p - ring[i];
		v.normalize();
		v0.normalize();
		v1.normalize();
		double angle[3];
		angle[0] = acos(-(v | v0)) * 180 / PI;
		angle[1] = acos((v | v1)) * 180 / PI;
		angle[2] = acos((v0 | v1)) * 180 / PI;
		for (int j = 0; j < 3; j++){
			if (minangle>angle[j]) minangle = angle[j];
			if (maxangle<angle[j]) maxangle = angle[j];
		}
	}
	return fabs(minangle - 0) + 0*fabs(maxangle - 60);
}

double computeLocalAngleQuality(MyTetMesh &mesh, OpenVolumeMesh::VertexHandle vh, Vec3d &cog){
	double minangle = 360;

	//遍历点周围的单元
	for (auto vc_it = mesh.vc_iter(vh); vc_it.valid(); ++vc_it){
		vector<Vec3d> vs;
		vs.push_back(cog);
		for (auto tv_it = mesh.cv_iter(*vc_it); tv_it.valid(); ++tv_it){
			if (*tv_it != vh)
				vs.push_back(mesh.vertex(*tv_it));
		}
		vector<Vec3d> normals;
		Vec3d v1;
		Vec3d v2;
		Vec3d v3;
		v1 = vs[1] - vs[0];
		v2 = vs[2] - vs[0];
		//v3 = vs[3] - vs[0];
		//if ((v3 | (v1%v2)) < 0)
		//	normals.push_back(v2%v1);
		//else
		normals.push_back(v1%v2);
		v1 = vs[2] - vs[0];
		v2 = vs[3] - vs[0];
		normals.push_back(v1%v2);
		v1 = vs[3] - vs[0];
		v2 = vs[1] - vs[0];
		normals.push_back(v1%v2);
		v1 = vs[3] - vs[1];
		v2 = vs[2] - vs[1];
		normals.push_back(v1%v2);
		for (int i = 0; i < normals.size(); ++i){
			normals[i].normalize();
		}
		vector<double> angles;
		for (int i = 0; i < normals.size() - 1; ++i){
			for (int j = i + 1; j < normals.size(); ++j){
				angles.push_back(acos((normals[i] | normals[j])) / PI*180.);
			}
		}
		for (int i = 0; i < angles.size(); ++i){
			if (angles[i] < minangle)
				minangle = angles[i];
		}
	}
	return minangle;
}


bool isImprovedLocally(MyMesh &mesh, MyMesh::VertexHandle vh, MyMesh::Point &cog){
	//记录每个节点周围节点
	vector<MyMesh::Point> ringpoints;
	MyMesh::Point origin = mesh.point(vh);
	for (auto vv_it = mesh.vv_cwbegin(vh); vv_it != mesh.vv_cwend(vh); ++vv_it)
	{
		ringpoints.push_back(mesh.point(*vv_it));
	}
	if (computeLocalAngleQuality(cog, ringpoints) >= computeLocalAngleQuality(origin, ringpoints))
		return true;
	else return false;
}
bool isImprovedLocally(MyTetMesh &mesh, OpenVolumeMesh::VertexHandle vh, Vec3d &cog){
	Vec3d origin = mesh.vertex(vh);
	if (computeLocalAngleQuality(mesh, vh, cog) >= computeLocalAngleQuality(mesh, vh, origin))
		return true;
	else return false;
}


//输出网格质量
void outMeshIdealElementQuality(vector<double> &qualities1)
{
	int num1, num2, num3, num4, num5, num6, num7, num8, num9, num10;
	num1 = num2 = num3 = num4 = num5= num6= num7= num8= num9= num10 = 0;
	double mean = 0.;
	double maxangle = 0.;
	double minangle = 1.;
	int num = qualities1.size();
	for (int i = 0; i < qualities1.size(); i++)
	{
		if (qualities1[i] <= 0.1)
		{
			num1++;
		}
		else if (qualities1[i] <= 0.2)
		{
			num2++;
		}
		else if (qualities1[i] <= 0.3)
		{
			num3++;
		}
		else if (qualities1[i] <= 0.4)
		{
			num4++;
		}
		else if (qualities1[i] <= 0.5)
		{
			num5++;
		}
		else if (qualities1[i] <= 0.6)
		{
			num6++;
		}
		else if (qualities1[i] <= 0.7)
		{
			num7++;
		}
		else if (qualities1[i] <= 0.8)
		{
			num8++;
		}
		else if (qualities1[i] <= 0.9)
		{
			num9++;
		}
		else
		{
			num10++;
		}
		mean += qualities1[i];
		if (qualities1[i] > maxangle)
		{
			maxangle = qualities1[i];
		}
		if (qualities1[i] < minangle)
		{
			minangle = qualities1[i];
		}
	}
	cout << "min_idel: " << minangle << " max_idel: " << maxangle << endl;
	cout << num << " " << mean / num << endl; 
	cout << num1 << " " << num2 << " " << num3 << " " << num4 << " " << num5 << " " << num6 << " " << num7 << " " << num8 << " " << num9 << " " << num10 << endl;
}
void outMeshAngleQuality(vector<double> &qualities2)
{
	int num1, num2, num3, num4, num5, num6, num7, num8, num9, num10, num11, num12, num13, num14, num15, num16, num17, num18;
	num1 = num2 = num3 = num4= num5= num6= num7= num8= num9= num10= num11= num12= num13= num14= num15= num16= num17= num18 = 0;
	double mean = 0.;
	double maxangle = 0.;
	double minangle = 180.;
	int num = qualities2.size();
	for (int i = 0; i < qualities2.size(); i++)
	{
		if (qualities2[i] <= 10.)
		{
			num1++;
		}
		else if (qualities2[i] <= 20.)
		{
			num2++;
		}
		else if (qualities2[i] <= 30.)
		{
			num3++;
		}
		else if (qualities2[i] <= 40.)
		{
			num4++;
		}
		else if (qualities2[i] <= 50.)
		{
			num5++;
		}
		else if (qualities2[i] <= 60.)
		{
			num6++;
		}
		else if (qualities2[i] <= 70.)
		{
			num7++;
		}
		else if (qualities2[i] <= 80.)
		{
			num8++;
		}
		else if (qualities2[i] <= 90.)
		{
			num9++;
		}
		else if (qualities2[i] <= 100.)
		{
			num10++;
		}
		else if (qualities2[i] <= 110.)
		{
			num11++;
		}
		else if (qualities2[i] <= 120.)
		{
			num12++;
		}
		else if (qualities2[i] <= 130.)
		{
			num13++;
		}
		else if (qualities2[i] <= 140.)
		{
			num14++;
		}
		else if (qualities2[i] <= 150.)
		{
			num15++;
		}
		else if (qualities2[i] <= 160.)
		{
			num16++;
		}
		else if (qualities2[i] <= 170.)
		{
			num17++;
		}
		else if (qualities2[i] > 170.)
		{
			num18++;
		}
		mean += qualities2[i];
		if (qualities2[i] > maxangle)
		{
			maxangle = qualities2[i];
		}
		if (qualities2[i] < minangle)
		{
			minangle = qualities2[i];
		}
	}
	num = qualities2.size();
	cout << "min_angle: " << minangle << " max_angle: " << maxangle << endl;
	cout << num << " " << mean / num << endl; 
	cout << num1 << " " << num2 << " " << num3 << " " << num4 << " " << num5 << " " << num6 << " " << num7 << " " << num8 << " " << num9 << " " << num10 << " " << num11 << " " << num12 << " " << num13 << " " << num14 << " " << num15 << " " << num16 << " " << num17 << " " << num18 << endl;
}
void outFileMeshAngleQuality(string fname, vector<double> &qualities2)
{
	fname.append(".txt");
	int num1, num2, num3, num4;
	num1 = num2 = num3 = num4 = 0;
	double mean = 0.;
	double maxangle = 0.;
	double minangle = 180.;
	int num = qualities2.size();
	for (int i = 0; i < qualities2.size(); i++)
	{
		if (qualities2[i] <= 10.)
		{
			num1++;
		}
		else if (qualities2[i] <= 20.)
		{
			num2++;
		}
		else if (qualities2[i] >= 165.)
		{
			num3++;
		}
		else if (qualities2[i] >= 150.)
		{
			num4++;
		}
		mean += qualities2[i];
		if (qualities2[i] > maxangle)
		{
			maxangle = qualities2[i];
		}
		if (qualities2[i] < minangle)
		{
			minangle = qualities2[i];
		}
	}
	ofstream fout;
	fout.open(fname,ios::app); 
	num = qualities2.size();
	fout << minangle << " " << maxangle << endl;
	//cout << num << " " << num1 << " " << num2 << " " << num3 << " " << num4 << " " << mean / num << endl;
}