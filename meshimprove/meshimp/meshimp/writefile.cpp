#include"writefile.h"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Handles.hh>
void writeVTKWithFixedBoundery(MyMesh &mesh, char* filename)
{
	string fname(filename); //[50];
	//cout << "output filename for vtkout\n" << endl;

	//cin>>fname;
	ofstream file(fname.c_str());
	file << "# vtk DataFile Version 3.0\n";
	file << "Mesquite Mesh\n";
	file << "ASCII\n";
	file << "DATASET UNSTRUCTURED_GRID\n";
	file << "POINTS " << mesh.n_vertices()<< " double\n";
	for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
	{
		file << mesh.point(*v_it)[0] << " " << mesh.point(*v_it)[1] << " " << mesh.point(*v_it)[2] << endl;
	}
	file << "CELLS " << mesh.n_faces() << " " << mesh.n_faces() * 4 << endl;
	for (auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
	{
		file << 3 ;
		for (auto fv_it = mesh.fv_begin(*f_it); fv_it != mesh.fv_end(*f_it); fv_it++)
		{
			file << " " << (*fv_it).idx();
		}
		file << endl;
	}
	file << "CELL_TYPES " << mesh.n_faces() << endl;
	for (int i = 0; i<mesh.n_faces(); i++)
		file << 5 << endl;
	file << "POINT_DATA " << mesh.n_vertices() << endl;
	file << "SCALARS fixed int\n";
	file << "LOOKUP_TABLE default\n";
	for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
	{
		if (mesh.is_boundary(*v_it))
			file << 1 << endl;
		else
			file << 0 << endl;
	}	
}

void readEleFile(MyTetMesh &mytet, string filename)
{
	string filename1(filename);
	filename1.append(".node");
	ifstream file1(filename1.c_str());
	int n, num;
	double a, b, c;
	file1 >> n>>a>>b>>c;
	vector<OpenVolumeMesh::VertexHandle> vhs;
	for (int i = 0; i < n; i++)
	{
		file1 >> num >> a >> b >> c;
		vhs.push_back(mytet.add_vertex(Vec3d(a, b, c)));
	}
	string filename2(filename);
	filename2.append(".ele");
	ifstream file2(filename2.c_str());
	int na, nb, nc, nd;
	file2 >> n >> na >> nb;
	for (int i = 0; i < n; i++)
	{
		file2 >> num >> na >> nb >> nc >> nd;
		mytet.add_cell(vhs[na-1],vhs[nb-1],vhs[nc-1],vhs[nd-1],1);
	}
}

void readQuadVTK(PolyMesh &mesh, char* filename)
{
	vector<PolyMesh::VertexHandle> vhs;
	string filena(filename);
	ifstream ist(filena.c_str());
	char buf[100];
	while (ist.good())
	{
		ist >> buf;

		int n = strlen(buf);
		for (int i = 0; i < n; i++)
			buf[i] = toupper(buf[i]);

		if (strcmp(buf, "POINTS") == 0)
		{
			int numv;
			ist >> numv;
			ist >> buf;
			for (int i = 0; i < numv; i++)
			{
				PolyMesh::Point point;
				ist >> point[0];
				ist >> point[1];
				ist >> point[2];
				vhs.push_back(mesh.add_vertex(point));
			}
		}

		if (strcmp(buf, "CELLS") == 0)
		{
			int a, b, c, d;
			int numc;
			ist >> numc;
			ist >> buf;
			for (int i = 0; i < numc; i++)
			{
				ist >> buf >> a >> b >> c >> d;
				vector<PolyMesh::VertexHandle> fhandle;
				fhandle.push_back(vhs[a]);
				fhandle.push_back(vhs[b]);
				fhandle.push_back(vhs[c]);
				fhandle.push_back(vhs[d]);
				mesh.add_face(fhandle);
				fhandle.clear();
			}
			break;
		}
	}
	cout << mesh.n_vertices() << endl;
	cout << mesh.n_faces() << endl;
}

void readTriVTK(MyMesh &mesh, char* filename)
{
	vector<MyMesh::VertexHandle> vhs;
	string filena(filename);
	ifstream ist(filena.c_str());
	char buf[100];
	while (ist.good())
	{
		ist >> buf;

		int n = strlen(buf);
		for (int i = 0; i < n; i++)
			buf[i] = toupper(buf[i]);

		if (strcmp(buf, "POINTS") == 0)
		{
			int numv;
			ist >> numv;
			ist >> buf;
			for (int i = 0; i < numv; i++)
			{
				PolyMesh::Point point;
				ist >> point[0];
				ist >> point[1];
				ist >> point[2];
				vhs.push_back(mesh.add_vertex(point));
			}
		}

		if (strcmp(buf, "CELLS") == 0)
		{
			int a, b, c;
			int numc;
			ist >> numc;
			ist >> buf;
			for (int i = 0; i < numc; i++)
			{
				ist >> buf >> a >> b >> c ;
				vector<MyMesh::VertexHandle> fhandle;
				fhandle.push_back(vhs[a]);
				fhandle.push_back(vhs[b]);
				fhandle.push_back(vhs[c]);
				mesh.add_face(fhandle);
				fhandle.clear();
			}
			break;
		}
	}
	cout << mesh.n_vertices() << endl;
	cout << mesh.n_faces() << endl;
}

void readTritxt(MyMesh &mesh, char* filename1, char* filename2)
{
	vector<MyMesh::VertexHandle> vhs;
	ifstream file1(filename1);
	int n, num;
	double a, b, c;
	file1 >> n;
	for (int i = 0; i < n; i++)
	{
		file1 >> a >> b;
		vhs.push_back(mesh.add_vertex(MyMesh::Point(a, b, 0.)));
	}
	ifstream file2(filename2);
	int na, nb, nc, nd;
	file2 >> n;
	for (int i = 0; i < n; i++)
	{
		file2 >> num >> na >> nb >> nc ;
		vector<MyMesh::VertexHandle> fhandle;
		fhandle.push_back(vhs[na-1]);
		fhandle.push_back(vhs[nb-1]);
		fhandle.push_back(vhs[nc-1]);
		mesh.add_face(fhandle);
		fhandle.clear();
	}
	cout << mesh.n_vertices() << endl;
	cout << mesh.n_faces() << endl;
}

void readliu(MyMesh &mesh, char* filename)
{
	vector<MyMesh::VertexHandle> vhs;
	ifstream file1(filename);
	int n, num;
	double a, b, c;
	int na, nb, nc;
	file1 >> n>>num;
	for (int i = 0; i < n; i++)
	{
		file1 >> a >> b >> c;
		vhs.push_back(mesh.add_vertex(MyMesh::Point(a, b, c)));
	}
	for (int i = 0; i < num; i++)
	{
		file1 >> na >> nb >> nc;
		vector<MyMesh::VertexHandle> fhandle;
		fhandle.push_back(vhs[na ]);
		fhandle.push_back(vhs[nb ]);
		fhandle.push_back(vhs[nc ]);
		mesh.add_face(fhandle);
		fhandle.clear();
	}
	cout << mesh.n_vertices() << endl;
	cout << mesh.n_faces() << endl;
}


void writeQuadVTK(PolyMesh &mesh, char* filename)
{
	string fname(filename); //[50];
	//cout << "output filename for vtkout\n" << endl;

	//cin>>fname;
	ofstream file(fname.c_str());
	file << "# vtk DataFile Version 3.0\n";
	file << "Mesquite Mesh\n";
	file << "ASCII\n";
	file << "DATASET UNSTRUCTURED_GRID\n";
	file << "POINTS " << mesh.n_vertices() << " double\n";
	for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
	{
		file << mesh.point(*v_it)[0] << " " << mesh.point(*v_it)[1] << " " << mesh.point(*v_it)[2] << endl;
	}
	file << "CELLS " << mesh.n_faces() << " " << mesh.n_faces() * 5 << endl;
	for (auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
	{
		file << 4;
		for (auto fv_it = mesh.fv_begin(*f_it); fv_it != mesh.fv_end(*f_it); fv_it++)
		{
			file << " " << (*fv_it).idx();
		}
		file << endl;
	}
	file << "CELL_TYPES " << mesh.n_faces() << endl;
	for (int i = 0; i<mesh.n_faces(); i++)
		file << 5 << endl;
	file << "POINT_DATA " << mesh.n_vertices() << endl;
	file << "SCALARS fixed int\n";
	file << "LOOKUP_TABLE default\n";
	for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
	{
		if (mesh.is_boundary(*v_it))
			file << 1 << endl;
		else
			file << 0 << endl;
	}
}

bool readTetVTK(MyTetMesh &mesh, string fname){
	fname.append(".vtk");
	ifstream fin(fname.c_str());
	string strtmp;
	for (; strtmp != "POINTS";) fin >> strtmp;
	int num_p;
	fin >> num_p >> strtmp;
	double x, y, z;
	vector<OpenVolumeMesh::VertexHandle> vhs;
	for (int i = 0; i < num_p; i++) {
		fin >> x >> y >> z;
		vhs.push_back(mesh.add_vertex(Vec3d(x, y, z)));
	}
	for (; strtmp != "CELLS";) fin >> strtmp;
	int num_c, vertex;
	fin >> num_c >> x;
	vector<OpenVolumeMesh::VertexHandle> vhs1;
	for (int i = 0; i < num_c; i++) {
		fin >> x;//4¸ö¶¥µã
		vhs1.clear();
		for (int j = 0; j < 4; j++) {
			fin >> vertex;
			vhs1.push_back(vhs[vertex]);
		}
		vector<OpenVolumeMesh::VertexHandle> vhs2;
		vhs2.push_back(vhs1[0]);
		vhs2.push_back(vhs1[1]);
		vhs2.push_back(vhs1[2]);
		vhs2.push_back(vhs1[3]);
		mesh.add_cell(vhs2);
	}
	fin.close();
	return true;
}