all:
# g++ -g -o fluid_sim main.cpp grid.cpp -lglut -lGL -lGLU
	g++ main.cpp fluidsolver.cpp -lGL -lGLU -lglut

clean:
	rm ./a.out