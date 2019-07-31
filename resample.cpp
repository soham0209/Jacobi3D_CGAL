#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <vector>
#include <Eigen/Dense>
#include <string>
#include <fstream>
#include <map>
typedef CGAL::Exact_predicates_exact_constructions_kernel         K;
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K>    Vb;
typedef CGAL::Delaunay_triangulation_cell_base_3<K>                 Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb>                Tds;
//Use the Fast_location tag. Default or Compact_location works too.
typedef CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location> Triangulation;
typedef Triangulation::Point	Point;
typedef Triangulation::Finite_cells_iterator Finite_cells_iterator;
typedef Triangulation::Finite_facets_iterator Finite_facets_iterator;
std::map<unsigned,double> f, g, h;
std::map<unsigned,Triangulation::Vertex_handle> v_info;
using std::cout; using std::endl;
int main(int argc, char* argv[])
{   std::string data = argv[1];
    std::vector< std::pair<Point,unsigned> > points;
    std::vector< std::pair<Point,unsigned> > added_points;
    std::string tri_file = data + "/" + data + "_triangulated.raw";
	std::string scalar_file = data + "/" + data + "_field.raw";
    std::string res_tri_file = data + "/" + data + "_resampled_triangulated.raw";
	std::string res_scalar_file = data + "/" + data + "_resampled_field.raw";
    std::string pers_file = data + "/" + data + "_resampled.txt";
    std::string vert_file_path =  data + "/" + data + "_resampled_vert.txt";
	std::ifstream tri(tri_file, std::ios::in);
	std::ifstream inStream(scalar_file, std::ios::binary);
    
	//std::ofstream out_scalar(scalar_file, std::ios::binary);
    Triangulation T;
	if (tri.good() && inStream.good()) {
		std::cout << "Triangulation already exists.Reading from it." << endl;
		tri >> T;
		std::cout << "Reading Done!" << std::endl;
		std::cout << "Triangulation consists " << T.number_of_vertices() << " Vertices and " << T.number_of_finite_cells() << " Cells." << endl;
		
		tri.close();
		cout << "Scalar values already exist.Reading from it." << endl;
		unsigned ind;
		double val;
		for (Triangulation::Finite_vertices_iterator vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); vit++) {
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
        unsigned id_max = T.number_of_vertices();
    
        //Triangulation T( points.begin(),points.end() );
        //CGAL_assertion( T.number_of_vertices() == 6 );
        // check that the info was correctly set.
        Triangulation::Finite_vertices_iterator vit;
        Triangulation::Finite_cells_iterator cit;
        for(cit = T.finite_cells_begin();cit!=T.finite_cells_end();cit++){
            //Triangulation::Cell_handle c = cit;
            bool greater_delta = false;
            double f_values = 0.0;
            double g_values = 0.0, h_values = 0.0;
            Point P = CGAL::centroid(T.tetrahedron(cit));
            for(int i=0; i<4; i++){
                Triangulation::Vertex_handle u = cit->vertex(i);
                if(!T.is_infinite(u)){
                    double f_val = f[u->info()];
                    double g_val = g[u->info()];
                    double h_val = h[u->info()];
                    
                    if( f_val > 10){
                        greater_delta = true;
                    }
                    f_values = f_values + f_val;
                    g_values = g_values + g_val;
                    h_values = h_values + h_val;
                }
            }
            //std::cout << v.point() << std::endl;
            if(greater_delta){
                Triangulation::Vertex v = cit->circumcenter();
                
                added_points.push_back( std::make_pair(P, id_max));
                f[id_max] = f_values/4.0;
                g[id_max] = f_values/4.0;
                h[id_max] = f_values/4.0;


                id_max ++;
            }
            //Triangulation::Vertex_handle vin = T.insert(v.point());
        }
        cout<<"Added points size "<< added_points.size() << endl;
        T.insert(added_points.begin(),added_points.end());
        std::cout << "Num Vertices " << T.number_of_vertices() << std::endl;
        std::cout << "Num Cells " << T.number_of_finite_cells() << std::endl;
        std::ofstream out_tri(res_tri_file, std::ios::out);
        //assert(g.size()==h.size()==f.size()==T.number_of_vertices());
        std::cout<<"F size " <<f.size()<<" gsize "<< g.size() << 
        " hsize "<< h.size() << endl; 
        out_tri << T;
        out_tri.close();
        std::ofstream ostream(res_scalar_file, std::ios::binary);
        std::ofstream pers_file_out(pers_file);
        std::ofstream vert_file(vert_file_path);
        pers_file_out << "3" << endl;
        pers_file_out << T.number_of_vertices() << endl;
        pers_file_out << "1" << endl;
        pers_file_out << "1" << endl;
		for (vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); vit++) {
			unsigned i = vit->info();
			ostream.write(reinterpret_cast<const char*>(&i), sizeof(i));
			ostream.write(reinterpret_cast<const char*>(&f[i]), sizeof(f[i]));
			ostream.write(reinterpret_cast<const char*>(&g[i]), sizeof(g[i]));
			ostream.write(reinterpret_cast<const char*>(&h[i]), sizeof(h[i]));
            v_info[i] = vit;
		}
		ostream.close();
        for(unsigned i =0;i<T.number_of_vertices();i++){
            assert(i == v_info[i]->info());
            pers_file_out << f[i] << endl;
            vert_file << v_info[i]->point() << endl;  
        }
        pers_file_out.close();
        vert_file.close();
		cout << "Saving Scalar Field Done" << endl;


        /*for (vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit){
            std::cout << "Point "<<vit->info()<<": "<<vit->point()<<std::endl;
        }
        std::cout << "OK" << std::endl;*/
        return 0;
    }
    return -1;
}
