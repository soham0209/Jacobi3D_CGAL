#include "Jacobi3d.h"
#include <Eigen/Dense>
#include <chrono>
using namespace std;
Jacobi3d::Jacobi3d(string filename) {
	vector<int> dirctn1, dirctn2;
	srand(time(0));
	int i = 0;
	while (i < 3) {
		int r = rand() % 9;
		if (r == 0)
			continue;
		dirctn1.push_back(r);
		i++;
	}
	dirctn2.push_back(dirctn1[2]);
	dirctn2.push_back(dirctn1[2]);
	dirctn2.push_back(-dirctn1[0] - dirctn1[1]);
	Eigen::Vector3d dir1(dirctn1[0], dirctn1[1], dirctn1[2]);
	Eigen::Vector3d dir2(dirctn2[0], dirctn2[1], dirctn2[2]);
	this->dir1 = dir1;
	this->dir2 = dir2;
	cout << "Dir1: " << dir1[0] << " " << dir1[1] << " " << dir1[2] << endl;
	cout << "Dir2: " << dir2[0] << " " << dir2[1] << " " << dir2[2] << endl;
	triangulate(filename);
}
void Jacobi3d::triangulate(string filename) {

	string tri_file = filename.substr(0, filename.find("_vert")) + "_triangulated.raw";
	string scalar_file = filename.substr(0, filename.find("_vert")) + "_field.raw";
	ifstream tri(tri_file, std::ios::in);
	ifstream inStream(scalar_file, std::ios::binary);
	
	if (tri.good() && inStream.good()) {
		cout << "Triangulation already exists.Reading from it." << endl;
		tri >> this->T;
		cout << "Reading Done!" << endl;
		cout << "Triangulation consists " << T.number_of_vertices() << " Vertices and " << T.number_of_finite_cells() << " Cells." << endl;
		
		tri.close();
		cout << "Scalar values already exist.Reading from it." << endl;
		unsigned ind;
		double val;
		for (auto vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); vit++) {
			inStream.read(reinterpret_cast<char*>(&ind), sizeof(ind));
			vit->info() = ind;
			inStream.read(reinterpret_cast<char*>(&val), sizeof(val));
			f[ind] = val;
			inStream.read(reinterpret_cast<char*>(&val), sizeof(val));
			g[ind] = val;
			inStream.read(reinterpret_cast<char*>(&val), sizeof(val));
			h[ind] = val;
		}
		inStream.close();
		this->read_from_saved = true;
	}
	else {
		ifstream inFile(filename);
		if (!inFile.is_open()) {
			cerr << "Input Triangulation not read" << endl;
		}
		else {
			list<std::pair<Point, unsigned>> L;
			unsigned int i = 0;
			double a, b, c;
			while (inFile >> a >> b >> c) {
				Point p(a, b, c);
				L.push_front(std::make_pair(p, i));
				Eigen::Vector3d v(a, b, c);
				g[i] = v.dot(dir1);
				h[i] = v.dot(dir2);
				i++;
			}
			cout << "Total vertices read " << i << endl;

			Triangulation K(L.begin(), L.end());
			cout << "Triangulation done !!!!" << endl;

			std::ofstream oFileT(tri_file, std::ios::out);
			// writing file output;
			oFileT << K;
			this->T = K;
		}
		inFile.close();
	}
}

void Jacobi3d::readperseus(string filename) {
	if (this->read_from_saved)
		return;
	cout << "Reading Scalar Field" << endl;
	string scalar_file = filename.substr(0, filename.find(".")) + "_field.raw";
	ifstream pers_file(filename);
	if (!pers_file.is_open()) {
		cerr << "Field not read" << endl;
	}
	else {
		string temp;
		pers_file >> temp;
		unsigned dim = stoi(temp);
		vector<unsigned> size_in_each_dim;
		unsigned total_ver = 1;
		for (unsigned i = 0; i < dim; i++) {
			int num_in_each;
			pers_file >> num_in_each;
			size_in_each_dim.push_back(num_in_each);
			total_ver = total_ver * num_in_each;
		}
		for (unsigned i = 0; i < total_ver; i++) {
			pers_file >> temp;
			double offset = (double)i / (double)total_ver;
			double val = stod(temp) + offset;
			f[i] = val;
			g[i] = g[i] + std::pow(offset, 2);
			h[i] = h[i] + std::pow(offset, 3);
			
		}
		pers_file.close();
		ofstream ostream(scalar_file, std::ios::binary);
		for (auto vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); vit++) {
			unsigned i = vit->info();
			ostream.write(reinterpret_cast<const char*>(&i), sizeof(i));
			ostream.write(reinterpret_cast<const char*>(&f[i]), sizeof(f[i]));
			ostream.write(reinterpret_cast<const char*>(&g[i]), sizeof(g[i]));
			ostream.write(reinterpret_cast<const char*>(&h[i]), sizeof(h[i]));
		}
		ostream.close();
		cout << "Reading and Saving Scalar Field Done" << endl;
	}
}

void Jacobi3d::getJacobiEdges() {
	for (auto fit = T.finite_facets_begin(); fit != T.finite_facets_end(); fit++) {
		vector<unsigned> link;
		vector<unsigned> face_ver;
		getLink(fit, link);
		if (link.size() < 2)
			continue;
		for (int i = 0; i < 4; i++) {
			if (i != fit->second) {
				Triangulation::Vertex_handle v = fit->first->vertex(i);
				face_ver.push_back(v->info());
			}
		}
		assert(link.size() == 2);
		assert(face_ver.size() == 3);
		if (isJacobi(face_ver, link)) {
			this->j_faces.push_back(face_ver);
		}
	}
}
void Jacobi3d::getLink(Finite_facets_iterator fit, vector<unsigned> &lk) {

	Triangulation::Cell_handle c = fit->first;
	Triangulation::Cell_handle c_nbr = c->neighbor(fit->second);
	if (!T.is_infinite(c_nbr) && !T.is_infinite(c)) {
		Triangulation::Vertex_handle u = c->vertex(fit->second);
		lk.push_back(u->info());
		for (int i = 0; i < 4; i++) {
			Triangulation::Vertex_handle v = c_nbr->vertex(i);
			if (!c->has_vertex(v)) {
				lk.push_back(v->info());
			}
		}
	}
}

bool Jacobi3d::isJacobi(vector<unsigned> face, vector<unsigned> lk) {
	bool is_lower_v1 = is_lower_link(face, lk[0]);
	bool is_lower_v2 = is_lower_link(face, lk[1]);
	return is_lower_v1 == is_lower_v2;
}

bool Jacobi3d::is_lower_link(vector<unsigned> face_ver, unsigned u) {
	Eigen::Vector3d v(f[u], g[u], h[u]);
	Eigen::Vector3d a(f[face_ver[0]], g[face_ver[0]], h[face_ver[0]]);
	Eigen::Vector3d b(f[face_ver[1]], g[face_ver[1]], h[face_ver[1]]);
	Eigen::Vector3d c(f[face_ver[2]], g[face_ver[2]], h[face_ver[2]]);
	Eigen::Vector3d ab = a - b;
	Eigen::Vector3d cb = c - b;
	Eigen::Vector3d va = v - a;
	return va.dot(ab.cross(cb)) < 0;
}

void Jacobi3d::WriteToFile(string filename) {
	ofstream fout(filename);
	cout << "Writing " << j_faces.size() << " faces ..." << endl;
	for (size_t i = 0; i < j_faces.size(); i++) {
		auto s = j_faces[i];
		for (auto it = s.begin(); it != s.end(); it++) {
			//cout << *it << " ";
			fout << *it << " ";
		}
		//cout << endl;
		fout << endl;
	}
	fout.close();
}
void Jacobi3d::debug() {
	for (auto vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); vit++) {
		unsigned i = vit->info();
		cout << i << " " << vit->point() <<" "<< f[i] << " " << g[i] << " " << h[i] << endl;
	}
}
