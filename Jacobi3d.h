#pragma once
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <map>
#include <set>
#include <ctime>
typedef CGAL::Exact_predicates_exact_constructions_kernel         K;
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K>    Vb;
typedef CGAL::Delaunay_triangulation_cell_base_3<K>                 Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb>                Tds;
//Use the Fast_location tag. Default or Compact_location works too.
typedef CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location> Triangulation;
typedef Triangulation::Point	Point;
typedef Triangulation::Finite_cells_iterator Finite_cells_iterator;
typedef Triangulation::Finite_facets_iterator Finite_facets_iterator;
#include <string>
using namespace std;
class Jacobi3d {
public:
	Jacobi3d(string filename);
	void triangulate(string filename);
	void readperseus(string filename);
	void getJacobiEdges();
	void WriteToFile(string filename);
	void debug();

private:	
Triangulation T;
void getLink(Finite_facets_iterator fit, vector<unsigned>& lk);
bool isJacobi(vector<unsigned> faces, vector<unsigned> lk);
bool is_lower_link(vector<unsigned> faces, unsigned u);
Eigen::Vector3d dir1, dir2;
public:
	map<unsigned,double> f, g, h;
	vector<vector<unsigned> > j_faces;
	bool read_from_saved = false;
	
};
