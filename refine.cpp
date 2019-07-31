#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Image_3.h>
#include <functional>
#include <boost/parameter.hpp>
#include <string>
#include <iterator>
typedef float Image_word_type;
// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Labeled_mesh_domain_3<K> Mesh_domain;
// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
// To avoid verbose function and named parameters call
using namespace CGAL::parameters;
int main(int argc, char*argv[])
{
  std::string data = argv[1];
  std::string im = data+"/"+data+".raw";
  const char* fname = im.c_str();
  const unsigned dim = (argc>2)?std::stoi(argv[2]):30;
  // Load image
  CGAL::Image_3 image;
  std::cout<<fname<<" "<<dim<<std::endl;
  if(!image.read_raw(fname,dim,dim,dim,1.0,1.0,1.0,0U,8,WK_FLOAT,SGN_SIGNED)){
    std::cerr << "Error: Cannot read file " <<  fname << std::endl;
    return EXIT_FAILURE;
  }
  Mesh_domain domain =
    Mesh_domain::create_gray_image_mesh_domain(image,value_outside=0.0);
  // Mesh criteria
  Mesh_criteria criteria(facet_angle=30, facet_size=6, facet_distance=2,
                         cell_radius_edge_ratio=3, cell_size=8);
  
  // Meshing
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);
  
  // Output
  std::ofstream vert_file(data+"/"+data+"_resampled_vert.txt");
  std::ofstream txt_file(data+"/"+data+"_resampled.txt");
  //c3t3.output_to_medit(medit_file);
  //c3t3.output_boundary_to_off(off_file);
  std::cout<<"V: "<<c3t3.number_of_vertices_in_complex()
    <<" F: "<<c3t3.number_of_facets()
    <<" C: "<<c3t3.number_of_cells()<<
    " C_in_cpx: "<<c3t3.number_of_cells_in_complex()<<std::endl;
    C3t3::Vertices_in_complex_iterator vit = c3t3.vertices_in_complex_begin();
    C3t3::Triangulation tr = c3t3.triangulation();
  int count = 1;
  int num_pt = std::distance(tr.points_begin(),tr.points_end());
  txt_file << "3" <<std::endl;
  txt_file <<  num_pt << std::endl;
  txt_file << "1" <<std::endl;
  txt_file << "1" <<std::endl;
  for(auto v = tr.points_begin();v!=tr.points_end();v++)
  {
    std::cout << "Pt "<<count<<": "<< v->point() 
    << " val: "<<image.value(v->x(),v->y(),v->z())<<std::endl;
    vert_file << v->point() << std::endl;
    txt_file << image.value(v->x(),v->y(),v->z()) << std::endl;
    count++;
  }
  //assert(num_pt == count);
  std::cout<<num_pt<<" "<<count<<std::endl;
  txt_file.close();
  vert_file.close();
  //std::cout<<sizeof(double)<<std::endl;
  /*
  std::cout<<"Size: "<<image.size()<<std::endl;
  for(size_t z = 0;z<dim;z++){
    for(size_t y = 0;y<dim;y++){
      for(size_t x = 0;x<dim;x++){
          auto val = image.value(x,y,z);
          if (val > 0){
            size_t ind  = (z * 4) + (y * 2) + x;
            std::cout << ind << " : " << val << std::endl; 
          }
        }
      }
    }*/
  return 0;
}