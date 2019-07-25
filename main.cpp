#include <iostream>
#include <string>
#include "Jacobi3d.h"

int main(int argc, char* argv[])
{
	/*
	// construction from a list of points :
	std::list<std::pair<Point,unsigned>> L;
	
	L.push_front(std::make_pair(Point(0, 0, 0),0));
	L.push_front(std::make_pair(Point(1, 0, 0),1));
	L.push_front(std::make_pair(Point(0, 1, 0),2));
	L.push_front(std::make_pair(Point(0, 1, 1),3));
	L.push_front(std::make_pair(Point(0, 1, -1),4));

	Triangulation T(L.begin(), L.end());
	
	Triangulation::Finite_vertices_iterator vit = T.finite_vertices_begin();
	for (; vit != T.finite_vertices_end(); vit++) {
		std::cout << vit->info() <<" "<< vit->point() << std::endl;
	}

	Finite_cells_iterator cit = T.finite_cells_begin();
	
	Triangulation::Finite_edges_iterator eit = T.finite_edges_begin();
	int count = 0;

	for (auto fit = T.finite_facets_begin(); fit != T.finite_facets_end(); fit++) {
		Triangulation::Cell_handle c = fit->first;
		Triangulation::Cell_handle c_nbr = c->neighbor(fit->second);
		if (!T.is_infinite(c_nbr) && !T.is_infinite(c)) {
			Triangulation::Vertex_handle u = c->vertex(fit->second);

			std::cout << "Lk Face " << count << ": (" << c->vertex(fit->second)->point() << ")"<< std::endl;
			std::cout << "Cell Info: " << std::endl;
			for (int i = 0; i < 4; i++) {
				std::cout << "Point " << i << ": " << c->vertex(i)->point() << std::endl;
			}
			std::cout << "Neighbor Info: " << std::endl;
			for (int i = 0; i < 4; i++) {
				Triangulation::Vertex_handle v = c_nbr->vertex(i);
				std::cout << "Point " << i << ": " << v->point() << std::endl;
				std::cout << c->has_vertex(v) << std::endl;
			}
			
			count++;
		}
	}*/
	if (argc < 2) {
		cout << "Argument missing" << endl;
		exit(1);
	}
	bool resamp = false;
	string data = argv[1];
	string resampled;
	if(argc > 2){
		resampled = argv[2];
		if(resampled == "-r")
			resamp = true;
		
	}
	string vert_file = data + "/" + data + "_vert.txt";
	string pers_file = data + "/" + data + ".txt";
	if(resamp){
		vert_file = data + "/" + data + "_resampled_vert.txt";
		pers_file = data + "/" + data + "_resampled.txt";

	}

	Jacobi3d J(vert_file);
	J.readperseus(pers_file);
	//J.debug();
	J.getJacobiEdges();
	J.WriteToFile(data + "/" + data + "_Jacobi.txt");
	return 0;
}
