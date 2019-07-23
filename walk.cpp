#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <map>
#include <set>
#include <ctime>
#include <string>
#include <algorithm>
#include <stack>
typedef CGAL::Exact_predicates_exact_constructions_kernel         K;
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K>    Vb;
typedef CGAL::Delaunay_triangulation_cell_base_3<K>                 Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb>                Tds;
//Use the Fast_location tag. Default or Compact_location works too.
typedef CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location> Triangulation;
typedef Triangulation::Point	Point;
typedef Triangulation::Finite_cells_iterator Finite_cells_iterator;
typedef Triangulation::Finite_facets_iterator Finite_facets_iterator;
using namespace std;

int readFaces(string dataFile, map<string,bool> &pruned_faces);
void traverse_triangulation(Triangulation T, map<string,bool> pruned_faces, vector<vector<int> > &faces_to_print);
string getface(Triangulation::Cell_handle u,Triangulation::Cell_handle v, vector<int> &faces);
int writefaces(string outfile, vector<vector<int> > faces_to_print);

int main(int argc, char* argv[]){
    string data = argv[1];
    string filename = data + "/" + data + "_vert.txt";
    string tri_file = filename.substr(0, filename.find("_vert")) + "_triangulated.raw";
	string scalar_file = filename.substr(0, filename.find("_vert")) + "_field.raw";
	ifstream tri(tri_file, std::ios::in);
	ifstream inStream(scalar_file, std::ios::binary);
	Triangulation T;
    map<unsigned, double> f, g, h;
    vector<vector<int> > faces_to_print;

	if (tri.good() && inStream.good()) {
		cout << "Triangulation already exists.Reading from it." << endl;
		tri >> T;
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
        map<string, bool> faces;
        int count = readFaces(filename, faces);
        cout << "Read " << count << " faces" << endl;
		if (count > 0) {
			traverse_triangulation(T, faces, faces_to_print);
			count = writefaces(data + "/" + data + "_pruned_walked.txt", faces_to_print);
			cout << "Wrote " << count << " faces to " << data + "/" + data + "_pruned_walked.txt" << endl;
		}
    }
    return 0;

}

int readFaces(string dataFile, map<string,bool> &pruned_faces) {
    string tri_file = dataFile.substr(0, dataFile.find("_vert")) + "_pruned_jacobi.off";
    ifstream pruned_file(tri_file);
    string temp;
    if (pruned_file.good()){
        pruned_file >> temp;
        if (temp!="OFF"){
            cout << "Bad OFF File. Header mismatch"<<endl;
            return -1;
        }
        pruned_file >> temp;
        int num_v = stoi(temp);
        pruned_file >> temp;
        int num_f = stoi(temp);
        pruned_file >> temp;
        for (int i =0;i<num_v;i++){
          int x,y,z;
          pruned_file >> temp;
          x = stoi(temp);
          pruned_file >> temp;
          y = stoi(temp);
          pruned_file >> temp;
          z = stoi(temp);  
        }
        for (int i =0;i<num_f;i++){
            pruned_file >> temp;
            int dim = stoi(temp);
            string faceid("");
            vector<unsigned> face_id;
            for(int j =0;j<dim;j++){
                pruned_file >> temp;
                face_id.push_back(stoi(temp));
            }
            sort(face_id.begin(),face_id.end());
            for(size_t i =0;i<face_id.size();i++){
                faceid = faceid + to_string(face_id[i]);
            }
            pruned_faces[faceid] = true;
        }
        return num_f;
    }
    cout << "OFF file not read."<< endl;
    return -1;
}

void traverse_triangulation(Triangulation T,map<string,bool> pruned_faces, vector<vector<int> > &faces_to_print){
    Finite_cells_iterator fit = T.finite_cells_begin();
    stack<Triangulation::Cell_handle> s;
    map<Triangulation::Cell_handle,bool> is_discovered;
    s.push(fit);
    while(!s.empty()){
        Triangulation::Cell_handle v = s.top();
        //cout<<s.size()<<endl;
        s.pop();
        if(!is_discovered[v] || is_discovered.count(v) == 0){
            is_discovered[v] = true;
            for(int i =0;i<4;i++){
                Triangulation::Cell_handle u = v->neighbor(i);
                if(!T.is_infinite(u)){
                    vector<int> common_face;
                    string f = getface(u, v, common_face);
                    if (pruned_faces.count(f)!= 0){
                        faces_to_print.push_back(common_face);   
                        continue;
                    }
                    else {
                        s.push(u);
                    }
                }
            }

        }

    }
    
}
string getface(Triangulation::Cell_handle u,Triangulation::Cell_handle v, vector<int> &faces){
    for (int i =0; i<4; i++){
        Triangulation::Vertex_handle p = u->vertex(i);
        if(v->has_vertex(p)){
            faces.push_back(p->info());
        }
    }
    sort(faces.begin(),faces.end());
    string faceid("");
    assert(faces.size() == 3);
    for(size_t i =0 ;i<faces.size();i++){
        faceid = faceid + to_string(faces[i]);
    }
    return faceid;
}
int writefaces(string outfile, vector<vector<int> > faces_to_print){
    ofstream ofile(outfile);
    if(ofile.is_open()){
        for(auto f:faces_to_print){
            for(auto v: f){
                ofile << v << " ";
            }
            ofile << endl;
        }
    ofile.close();
    //cout << "Wrote " << faces_to_print.size() << " to " << outfile << endl; 
    return (int)faces_to_print.size();
    }
    return -1;
    
}