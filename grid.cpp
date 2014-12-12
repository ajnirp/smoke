#include "grid.h"

Grid_Cell::Grid_Cell(int a_x, int a_y)
{
	x = a_x;
	y = a_y;
	vel_x = vel_y =0;

	//present = false;

	pressure = 0;

	for(int i=0; i< 5 ;i++){
		boundary[i] = false;
	}
}

Grid_Cell::~Grid_Cell(){

}

bool Grid_Cell::particle_present()
{
	return !(contained_particles.empty());
}

void Grid_Cell::clear_vals(){
	vel_x = 0;
	vel_y = 0;

	temp_count_x = 0;
	temp_count_y = 0;

	/*for(int i=0;i<5;i++){
		boundary[i] = false;
	}*/
}

Particle::Particle(float a_x, float a_y, int a_id)
{
	x_pos = a_x;
	y_pos = a_y;

	vel_x = vel_y = 0;

	id = a_id;
}

Particle::~Particle(){

}

Grid::Grid(float a_t, float a_cell_size)
{
	cell_size = a_cell_size;
	time_step = a_t;
}

Grid::~Grid(){}

void Grid::initialise_grid(int a_w, int a_h)
{
	w = a_w;
	h = a_h;

	for(int i = 0; i< w+1 ; i++){
		vector < Grid_Cell* > temp;
		temp.clear();
		for(int j = 0 ; j< h+1 ; j++)
		{
			temp.push_back(new Grid_Cell(i,j));
		}
		cells.push_back(temp);
	}

    //setting boundaries
	for(int i=0;i<w+1;i++){
		cells[i][0]->boundary[4] = true;
		cells[i][h-1]->boundary[3] = true;
	}

	for(int j=0;j<h+1;j++){
		cells[0][j]->boundary[1] = true;
		cells[w-1][j]->boundary[2] = true;
	}
	//cout<<cells[0][0]->boundary[1]<<endl;

	//Initiaze previous cells for FLIP
	if(model == FLIP){
		for(int i = 0; i< w+1 ; i++){
			vector < Grid_Cell* > temp;
			temp.clear();
			for(int j = 0 ; j< h+1 ; j++)
			{
				temp.push_back(new Grid_Cell(i,j));
			}
			prev_cells.push_back(temp);
		}		
	}
}

void Grid::initialise_fluid(){

	int count = 0;

	for(int i=0;i<w/2; i++){
		for(int j = 0 ; j< h/2 ; j++)
		{
			Grid_Cell* curr_grid = cells[i][j];

			//cout<<w/2<<" y3ah"<<endl;
			// Initialising 4 particles in each grid cell for now;
			count++;
			Particle * new_p = new Particle(i+0.25,j+0.25,count);
			new_p->in_grid = curr_grid;
			curr_grid->contained_particles.push_back(new_p);
			particles.push_back(new_p);

			count++;
			new_p = new Particle(i+0.75,j+0.25,count);
			new_p->in_grid = curr_grid;
			curr_grid->contained_particles.push_back(new_p);			
			particles.push_back(new_p);

			count++;
			new_p = new Particle(i+0.25,j+0.75,count);
			new_p->in_grid = curr_grid;
			curr_grid->contained_particles.push_back(new_p);
			particles.push_back(new_p);

			count++;			
			new_p = new Particle(i+0.75,j+0.75,count);
			new_p->in_grid = curr_grid;
			curr_grid->contained_particles.push_back(new_p);
			particles.push_back(new_p);			

		}
	}


}


void Grid::advect_particles(){
	for(list<Particle* >::iterator it = particles.begin() ; it != particles.end() ; it++)
	{
		float old_xpos = (*it)->x_pos;
		float old_ypos = (*it)->y_pos;
		(*it)->x_pos = old_xpos + time_step * (*it)->vel_x;
		(*it)->y_pos = old_ypos + time_step * (*it)->vel_y;


		// boundary conditions for grid boundary
		// should change for general boundaries
		if((*it)->x_pos < 0) {
			(*it)->x_pos = -((*it)->x_pos);
			(*it)->vel_x = -((*it)->vel_x );
		}	
		if((*it)->y_pos < 0) {
			(*it)->y_pos = -((*it)->y_pos);
			(*it)->vel_y = -((*it)->vel_y );
		}	
		if((*it)->y_pos > h) {
			(*it)->y_pos = 2*h-((*it)->y_pos);
			(*it)->vel_y = -((*it)->vel_y );
		}			
		if((*it)->x_pos > w) {
			(*it)->x_pos = 2*w-((*it)->x_pos);
			(*it)->vel_x = -((*it)->vel_x );
		}	

		if((*it)->id == 100 ){
			cout<<(*it)->x_pos<<" "<<(*it)->y_pos<<" "<<(*it)->vel_x<<" "<<(*it)->vel_y<<"\n";
		}


		// Updating if the particle changes gridcell
		if ( (int) old_xpos  != (int) ((*it)->x_pos) || (int) old_ypos  != (int) ((*it)->y_pos) ){
			//cout<<"1\n";

			Grid_Cell* old_grid = (*it)->in_grid;
			Grid_Cell* new_grid = cells[(int) ((*it)->x_pos)][(int) ((*it)->y_pos) ];

			old_grid->contained_particles.remove(*it);
			new_grid->contained_particles.push_back(*it);

			(*it)->in_grid = new_grid;
		}
	}
}


void Grid::project_particles_to_grid(){

	for(int i=0;i<w+1;i++){
		for (int j = 0; j < h+1; j++)
		{
			Grid_Cell* curr = cells[i][j];
			curr->clear_vals();
		}
	}


	for(list<Particle* >::iterator it = particles.begin() ; it != particles.end() ; it++){

		Grid_Cell* curr_grid = (*it)->in_grid;

		int x_val = curr_grid->x;
		int y_val = curr_grid->y;

		float c_x_vel = (*it)->vel_x;
		float c_y_vel = (*it)->vel_y;

		float hx = (*it)->x_pos - (int)((*it)->x_pos);
		float hy = (*it)->y_pos - (int)((*it)->y_pos);

		// Bilenear interploation on MAC grid
		// At a grid cell, bottom and left faces store the corresponding velocities

		cells[x_val][y_val]->vel_x += (1-hx)*c_x_vel;
		cells[x_val][y_val]->vel_y += (1-hy)*c_y_vel;

		// keeping track of the count of the particles
		cells[x_val][y_val]->temp_count_x+=(1-hx);
		cells[x_val][y_val]->temp_count_y+=(1-hy);


		cells[x_val+1][y_val]->vel_x += (hx)*c_x_vel;

		cells[x_val+1][y_val]->temp_count_x += (hx);


		cells[x_val][y_val+1]->vel_y += (hy)*c_y_vel;

		cells[x_val][y_val+1]->temp_count_y += (hy);

		//cout<<c_y_vel<<endl;

	}

	for(int i=0;i<w;i++){
		for (int j = 0; j < h; j++)
		{
			//weighted average based on the number of particles present
			//if 0, particles, => vel = 0.
			if(cells[i][j]->temp_count_x != 0){
				cells[i][j]->vel_x = (cells[i][j]->vel_x / (float) cells[i][j]->temp_count_x);
			}
			else{
				cells[i][j]->vel_x = 0;
			}

			if(cells[i][j]->temp_count_y != 0){
				cells[i][j]->vel_y = (cells[i][j]->vel_y / (float) cells[i][j]->temp_count_y);
			}	
			else{
				cells[i][j]->vel_y = 0;
			}
		}
	}
}

//Currently, only gravity
void Grid::buoyant_force()
{
	for(int i = 0; i < w ; i++)
	{
		for (int j = 0; j < h; j++)
		{
			Grid_Cell* curr = cells[i][j];
			if (curr->particle_present) {
				// buoyancy is an upward force, so change y velocity vector
				curr->vel_y += time_step * buoyancy_alpha * (curr->temperature - temperature_air);
			}
		}
	}
}

int Grid::convert_grid_to_array_index(int i,int j){
	return (j)*w + i;
	//from 0 to w*h -1
}

void Grid::projection_and_pressure_solve(){

	SparseMatrix<double> A(h*w,h*w);
	std::vector<T> tripletList;
	VectorXd b(h*w);
	VectorXd x(h*w);
	float sq_cell_size= cell_size*cell_size;
	float c_c; // current cell presure coeff
	float t_c,b_c,l_c,r_c; //pressures coeffs in other cells

	float rhs;
	for(int i=0;i<w;i++){
		for(int j=0;j<h;j++){

			int curr = convert_grid_to_array_index(i,j);
			int top_index = convert_grid_to_array_index(i,j+1);
			int btm_index = convert_grid_to_array_index(i,j-1);
			int left_index = convert_grid_to_array_index(i-1,j);
			int right_index = convert_grid_to_array_index(i+1,j);
			
			
			//TODO: need to change the boundary conditions to include solid blocks in the middle

			Grid_Cell* curr_cell = cells[i][j];
			//for air and solid cells, no equations
			if(!(curr_cell->particle_present()) ||  curr_cell->boundary[0] ){
				b(curr) = 0;
				x(curr) = 0;
				continue;
			}

			c_c = t_c = b_c = l_c = r_c = rhs = 0;

			//cout<<curr_cell->boundary[1]<<endl;
			if(curr_cell->boundary[1]){
				rhs -= cells[i][j]->vel_x;
				c_c -= 1;
				l_c += 1;
			}
			else if(!(cells[i-1][j]->particle_present())){
				l_c += 1;
			}

			if(curr_cell->boundary[2]){
				rhs += cells[i+1][j]->vel_x;
				c_c -= 1;
				r_c += 1;
			}
			else if(!(cells[i+1][j]->particle_present())){
				r_c += 1;
			}

			if(curr_cell->boundary[3]){
				rhs += cells[i][j+1]->vel_y;
				c_c -= 1;
				t_c += 1;
			}
			else if(!(cells[i][j+1]->particle_present())){
				t_c += 1;
			}

			if(curr_cell->boundary[4]){
				rhs -= cells[i][j]->vel_y;
				c_c -= 1;
				b_c += 1;
			}
			else if(!(cells[i][j-1]->particle_present())){
				b_c += 1;
			}									
		
			c_c +=4;
			l_c--;r_c--;b_c--;t_c--;

			// adding to the triplet list
			if(c_c != 0){
				tripletList.push_back(T(curr,curr,c_c/sq_cell_size));
				//cout<<"C:"<<curr<<endl;
			}
			if(t_c != 0){
				tripletList.push_back(T(curr,top_index,t_c/sq_cell_size));
				//cout<<"t:"<<top_index<<endl;
			}
			if(b_c != 0){
				tripletList.push_back(T(curr,btm_index,b_c/sq_cell_size));
				//cout<<"b:"<<btm_index<<endl;
			}
			if(l_c != 0){
				tripletList.push_back(T(curr,left_index,l_c/sq_cell_size));
				//cout<<"l:"<<left_index<<endl;
			}
			if(r_c != 0){
				tripletList.push_back(T(curr,right_index,r_c/sq_cell_size));
				//cout<<"r:"<<right_index<<endl;
			}			

			// adding RHS to the vector 'b'
			rhs += -((cells[i+1][j]->vel_x - cells[i][j]->vel_x) + (cells[i][j+1]->vel_y - cells[i][j]->vel_y));

			rhs = rhs / cell_size;

			b(curr) = rhs;
		}
		
	}

	//	generating sparse matrix and solve for pressure
	A.setFromTriplets(tripletList.begin(), tripletList.end());

	MatrixXd dmat;

	dmat = MatrixXd(A);

//	cout<<dmat<<endl;

//	cout<<b<<endl;

	// solve Ax = b
	ConjugateGradient<SparseMatrix<double> > solver;
	solver.compute(A);
	if(solver.info()!=Success) {
		cout<<"decomposition during pressure solve CG failed";
		return;
	}
	x = solver.solve(b);
	if(solver.info()!=Success) {
		cout<<"solving for pressure by CG failed";
		sleep(10);
		return;
	}
	// solve for another right hand side:
	//x1 = solver.solve(b1);
	for (int i=0;i<h*w;i++){
		//cout<<x(i)<<"\n";
	}
	//sleep(1);
	cout<<"done\n";
	//sleep(10);


	//set pressures
	for(int i=0;i<w;i++){
		for(int j=0;j<h;j++){

			Grid_Cell* curr_cell = cells[i][j];
			if(!(curr_cell->particle_present()) ||  curr_cell->boundary[0] ){
				curr_cell->pressure = 0;
				
			}
			else{
				int index = convert_grid_to_array_index(i,j);
				curr_cell->pressure = x(index);
			}

		}	

	}

	// print info
	for(int j=h-1;j>=0;j--){
		for(int i=0;i<w;i++){
			//cout<<cells[i][j]->pressure<<" ";
		}
		//cout <<":\n";
	}		

	//	update velocities with new pressures
	for(int i=0;i<w;i++){
		for(int j=0;j<h;j++){	
			Grid_Cell* curr_cell = cells[i][j];

			if(curr_cell->boundary[1]){
				curr_cell->vel_x = 0;
			}
			else{
				curr_cell->vel_x -= ((curr_cell->pressure - cells[i-1][j]->pressure)/cell_size);
			}

			if(curr_cell->boundary[4]){
				curr_cell->vel_y = 0;
			}
			else{
				curr_cell->vel_y -= ((curr_cell->pressure - cells[i][j-1]->pressure)/cell_size);
			}

			// other cases TODO ?

		}
	}


}

void Grid::project_grid_to_particles_PIC(){

	for(list<Particle* >::iterator it = particles.begin() ; it != particles.end() ; it++){

		Grid_Cell* curr_grid = (*it)->in_grid;

		int x_val = curr_grid->x;
		int y_val = curr_grid->y;


		float hx = (*it)->x_pos - (int)((*it)->x_pos);
		float hy = (*it)->y_pos - (int)((*it)->y_pos);

		// Bilenear interploation on to particles

		(*it)->vel_x =  cells[x_val][y_val]->vel_x*(1-hx) + cells[x_val+1][y_val]->vel_x*(hx);
		(*it)->vel_y =  cells[x_val][y_val]->vel_y*(1-hy) + cells[x_val][y_val+1]->vel_y*(hy);

		if((*it)->id == 16 ){
			cout<<"proj g->p: "<<(*it)->x_pos<<" "<<(*it)->y_pos<<" "<<(*it)->vel_x<<" "<<(*it)->vel_y<<"\n";
		}

	}

}


void copy_grid_cell_velocities(Grid_Cell* src, Grid_Cell* dest){
	dest->vel_x = src->vel_x;
	dest->vel_y = src->vel_y;
}

void Grid::save_grid(){

	for(int i=0;i<w+1;i++){
		for (int j = 0; j < h+1; j++)
		{
			copy_grid_cell_velocities(cells[i][j],prev_cells[i][j]);
		}
	}

}

void Grid::project_grid_to_particles_FLIP(){

	for(list<Particle* >::iterator it = particles.begin() ; it != particles.end() ; it++){

		Grid_Cell* curr_grid = (*it)->in_grid;

		int x_val = curr_grid->x;
		int y_val = curr_grid->y;


		float hx = (*it)->x_pos - (int)((*it)->x_pos);
		float hy = (*it)->y_pos - (int)((*it)->y_pos);

		// Bilenear interploation on to particles but only interpolate and add the difference in FLIP

		(*it)->vel_x +=  (cells[x_val][y_val]->vel_x - prev_cells[x_val][y_val]->vel_x)*(1-hx) + (cells[x_val+1][y_val]->vel_x - prev_cells[x_val+1][y_val]->vel_x)*(hx);
		(*it)->vel_y +=  (cells[x_val][y_val]->vel_y - prev_cells[x_val][y_val]->vel_y)*(1-hy) + (cells[x_val][y_val+1]->vel_y - prev_cells[x_val][y_val+1]->vel_y)*(hy);

		if((*it)->id == 16 ){
			cout<<"proj g->p: "<<(*it)->x_pos<<" "<<(*it)->y_pos<<" "<<(*it)->vel_x<<" "<<(*it)->vel_y<<"\n";
		}

	}

}


void Grid::Exec_Time_Step(){
	advect_particles();
	project_particles_to_grid();

	if(model == FLIP){
		save_grid();
	}

	add_force();

	projection_and_pressure_solve();

	if(model == PIC){
		project_grid_to_particles_PIC();
	}
	else if(model == FLIP){
		project_grid_to_particles_FLIP();  //TODO
	}
	//sleep(3);
}


void Grid::Draw(){

	//float cell_size = (float) 10 / (float) h ;
	//cout<<cell_size<<"\n";

	for(int i=0;i<w;i++){
		for (int j = 0; j < h; j++)
		{
			if(cells[i][j]->particle_present()){
				glColor3f(0.0,0.0,1.0);
				glBegin(GL_QUADS);
					glVertex3f(i*cell_size,j*cell_size,0);
					glVertex3f((i+1)*cell_size,j*cell_size,0);
					glVertex3f((i+1)*cell_size,(j+1)*cell_size,0);
					glVertex3f(i*cell_size,(j+1)*cell_size,0);
				glEnd();
			}
			else{
				glColor3f(1.0,1.0,1.0);
				glBegin(GL_QUADS);
					glVertex3f(i*cell_size,j*cell_size,0);
					glVertex3f((i+1)*cell_size,j*cell_size,0);
					glVertex3f((i+1)*cell_size,(j+1)*cell_size,0);
					glVertex3f(i*cell_size,(j+1)*cell_size,0);
				glEnd();
			}

		}
	}	
}

