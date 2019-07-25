CXX = g++
LDLIBS  = -lCGAL -lgmp -lmpfr
CXXFLAGS = -std=c++11
all : Jacobi Walk clean-obj

Jacobi:Jacobi3d.o main.o
	$(CXX) $(CXXFLAGS) -o Jacobi $(LDLIBS) Jacobi3d.o main.o
Walk:CGAL_Walk.o
	$(CXX) $(CXXFLAGS) -o Walk $(LDLIBS) CGAL_Walk.o
clean-obj:
	$(RM) *.o *~
clean:
	$(RM) Jacobi Walk *.o *~