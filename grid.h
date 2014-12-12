#include <GL/glut.h>
#include <GL/gl.h>
#include <vector>
#include <list>
#include <iostream>

//eigen lib for conjugate gradient
#include "Eigen/IterativeLinearSolvers"

using namespace std;
using namespace Eigen;

typedef Eigen::Triplet<double> T;

class Grid_Cell;

class Particle
{
public:
	Particle(float a_x , float a_y, int a_id );
	~Particle();

	int id;

	float x_pos;
	float y_pos;

	Grid_Cell* in_grid;

	float vel_x;
	float vel_y;
};

class Grid_Cell{

public:
	Grid_Cell(int a_x, int a_y);
	~Grid_Cell();

	int x;	
	int y;

	float vel_x; // for the cell (i,j) the velocity stored is u(i+0.5,j)
	float vel_y; // for the cell (i,j) the velocity stored is u(i,j+0.5)
	//float vel_z;

	float temp_count_x;
	float temp_count_y;

	list<Particle*> contained_particles;

	bool particle_present();  // Indicator of whether a particle is present or not
	void clear_vals();

	// defined at the voxel center
	float temperature;
	float density;
	float implicit_surface;
	float pressure;

	bool boundary[5]; // 0th index = true => empty else solid ,rest for boundary  1 =  left , 2 = right , 3 = top , 4 = bottom
};

class Grid
{
public:
	//takes timestep and cellsize as arguments
	Grid(float a_t, float a_cell_size);
	~Grid();
	int w, h;
	float cell_size;
	float time_step;

	// pg 5, equation 14
	// ambient temperature of surrounding air, assumed constant
	float temperature_air;
	float buoyancy_alpha;

	vector< vector< Grid_Cell* > > cells;
	list <Particle* > particles;
	vector< vector< Grid_Cell* > > prev_cells;	// used for flip

	void initialise_grid(int a_w, int a_h);
	void initialise_fluid();

	void advect_particles();

	void project_particles_to_grid();

	void save_grid(); 	// for FLIP

	// buoyant force, pg 5, equation 14
	float buoyant_force();

	void check_boundary_conditions();
	void projection_and_pressure_solve();

	void project_grid_to_particles_PIC();
	void project_grid_to_particles_FLIP();

	void Exec_Time_Step();
	void Draw();

	//helpers
	int convert_grid_to_array_index(int i,int j);
};
