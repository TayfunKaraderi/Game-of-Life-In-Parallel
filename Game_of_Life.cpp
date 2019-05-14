#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <time.h>
#include <sstream>
#include <fstream>
#include <vector>

using namespace std;

//initialization
int N = 50; // Each process will have an N by N grid -- Note this is not the total size of the grid
int max_steps = 10; // number of iterations for game of life
bool periodic = false;
int test = 1; // used in the non periodic case only -- this is to check if a process has neighbour - if not, will send an array of zeros as ghost points.


//------------------- Game of Life Serial -----------------------
// will be used in the parallel region
int imax = N + 2, jmax = N + 2; // grid dimentions includes ghost points passed from other processes
vector< vector < bool > > grid, new_grid;
// create a lattice of size NxN for each process

int num_neighbours(int ii, int jj)
{
	int ix, jx;
	int cnt = 0;
	for (int i = -1; i <= 1; i++)
		for (int j = -1; j <= 1; j++)
			if (i != 0 || j != 0)
			{
				ix = (i + ii + imax) % imax;
				jx = (j + jj + jmax) % jmax;
				if (grid[ix][jx]) cnt++;
			}
	return cnt;
}

void grid_to_file(int it, int process_row, int process_column)
{
	stringstream fname;
	fstream f1;

	fname << "output" << "_" << "iteration_" << it << "_processrow_" << process_row << "_processcolumn_" << process_column << ".txt";
	f1.open(fname.str().c_str(), ios_base::out);
	for (int i = 1; i < imax-1; i++) //dont output boundaries received from neighbour processes
	{
		for (int j = 1; j < jmax-1; j++)
			f1 << grid[i][j] << "\t";
		f1 << endl;
	}
	f1.close();
}

void do_iteration(void)
{
	for (int i = 0; i < imax; i++)
		for (int j = 0; j < jmax; j++)
		{
			new_grid[i][j] = grid[i][j];
			int num_n = num_neighbours(i, j);
			if (grid[i][j])
			{
				if (num_n != 2 && num_n != 3)
					new_grid[i][j] = false;
			}
			else if (num_n == 3) new_grid[i][j] = true;
		}
	grid.swap(new_grid);
}
//------------------- End Game of Life Serial --------------------


//----------------------------------------------------------------
int id, p, tag_num = 1;

void find_dimensions(int p, int &rows, int &columns)		//A bit brute force - this can definitely be made more efficient!
{
	int min_gap = p;
	for (int i = 1; i <= p / 2; i++)
	{
		if (p%i == 0)
		{
			int gap = abs(p / i - i);

			if (gap < min_gap)
			{
				min_gap = gap;
				rows = i;
				columns = p / i;
			}
		}
	}

	if (id == 0)
		cout << "Divide " << p << " into " << rows << " by " << columns << " grid" << endl;
}

int rows, columns;
int id_row, id_column;

void id_to_index(int id, int &id_row, int &id_column)
{
	id_column = id % columns;
	id_row = id / columns;
}

int id_from_index(int id_row, int id_column)
{
	
	if (periodic == true) // find periodic neighbours id
	{
		return ((id_row + rows) % rows) * columns + (id_column + columns) % columns;
	}
	if (periodic == false) {
		if (id_row >= rows || id_row < 0)
		{
			test = -1;
		}
		if (id_column >= columns || id_column < 0)
		{
			test = -1;
		}

		//return id_row * columns + id_column;
		return ((id_row + rows) % rows) * columns + (id_column + columns) % columns;

	}

}


int main(int argc, char *argv[])
{

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	srand(time(NULL) + id * 10);

	find_dimensions(p, rows, columns);

	id_to_index(id, id_row, id_column);

	int cnt = 0;

	MPI_Request* request = new MPI_Request[8];

	grid.resize(imax, vector<bool>(jmax));
	new_grid.resize(imax, vector<bool>(jmax));

	// fill the lattice randomly
	srand(time(NULL)+1000*id);// seed by time
	int b; // random value of GoL lattice points
	for (int i = 1; i < N+1; i++) {
		for (int j = 1; j < N+1; j++) {
			//srand(rand() + 1000*i + 450*j);
			int a = rand() % 10;
			if (a < 5) { b = 0; }
			if (a >= 5) { b = 1; }
			grid[i][j] = b; //lattice values filled here
		}
	}

	// do n iterations
	for (int n = 0; n < max_steps; n++)
	{
		// -----values to send and receive-----
		//slice edge columns
		bool * left_edge = new bool[N]; 	bool * right_edge_recv = new bool[N];
		bool * right_edge = new bool[N];    bool * left_edge_recv = new bool[N];


		for (int i = 0; i < N; i++) { left_edge[i] = grid[i + 1][1]; }
		for (int i = 0; i < N; i++) { right_edge[i] = grid[i + 1][N]; }

		//slice edge rows
		bool * top_edge = new bool[N];     	 bool * bottom_edge_recv = new bool[N];
		bool * bottom_edge = new bool[N];    bool * top_edge_recv = new bool[N];
		for (int i = 0; i < N; i++) { top_edge[i] = grid[1][i + 1]; } //note: 0th row and column is the ghost points received from neighbours.
		for (int i = 0; i < N; i++) { bottom_edge[i] = grid[N][i + 1]; } //note: (N+1)th row and colunmns are the the ghost points received from neighbours.

		//non periodic and no neighbours - then, send an array of zeros as the ghost points
		bool *zeros = new bool[N]; for (int i = 0; i < N; i++) { zeros[i] = 0; }

		//edgevalues
		bool * left_top = new bool[1]; 	bool * left_bottom = new bool[1]; bool * right_top = new bool[1];  bool * right_bottom = new bool[1];
		left_top[0] = grid[1][1];	   bool *left_top_recv = new bool[1];
		left_bottom[0] = grid[N][1];  bool *left_bottom_recv = new bool[1];
		right_top[0] = grid[1][N];    bool *right_top_recv = new bool[1];
		right_bottom[0] = grid[N][N]; bool *right_bottom_recv = new bool[1];

		//non periodic and no neighbours - then, send zero as the ghost point
		bool *zero = new bool[1]; zero[0] = 0;


		//----------------------


		//sends and receives
		//----------
		//1. send left
		id_to_index(id, id_row, id_column);
		int idsend = id_from_index(id_row, id_column - 1); // if (idsend == -1) {idsend = id_from_index(id_row, columns - 1);}
		if (periodic == true || test == 1) { MPI_Isend(left_edge, N, MPI_BYTE, idsend, tag_num = 1, MPI_COMM_WORLD, &request[0]);}
		else { MPI_Isend(zeros, N, MPI_BYTE, idsend, tag_num = 1, MPI_COMM_WORLD, &request[0]); } // if not periodic and no neighbour for the process send zeros
		//receive from right
		id_to_index(id, id_row, id_column);
		int idrecv = id_from_index(id_row, id_column + 1); //if (idrecv == -1) {idrecv = id_from_index(id_row, 0);
		MPI_Irecv(right_edge_recv, N, MPI_BYTE, idrecv, tag_num = 1, MPI_COMM_WORLD, &request[0]);

		//2. send right
		id_to_index(id, id_row, id_column);
		idsend = id_from_index(id_row, id_column + 1); // if (idsend == -1) {idsend = id_from_index(id_row, columns - 1);}
		if (periodic == true || test == 1) { MPI_Isend(right_edge, N, MPI_BYTE, idsend, tag_num = 2, MPI_COMM_WORLD, &request[1]); }
		else { MPI_Isend(zeros, N, MPI_BYTE, idsend, tag_num = 2, MPI_COMM_WORLD, &request[1]); }
		//receive from left
		id_to_index(id, id_row, id_column);
		idrecv = id_from_index(id_row, id_column - 1); //if (idrecv == -1) {idrecv = id_from_index(id_row, 0);
		MPI_Irecv(left_edge_recv, N, MPI_BYTE, idrecv, tag_num = 2, MPI_COMM_WORLD, &request[1]);

		//3. send top
		id_to_index(id, id_row, id_column);
		idsend = id_from_index(id_row - 1, id_column); // if (idsend == -1) {idsend = id_from_index(id_row, columns - 1);}
		if (periodic == true || test == 1) { MPI_Isend(top_edge, N, MPI_BYTE, idsend, tag_num = 3, MPI_COMM_WORLD, &request[2]); }
		else { MPI_Isend(zeros, N, MPI_BYTE, idsend, tag_num = 3, MPI_COMM_WORLD, &request[2]); }
		//receive from bottom
		id_to_index(id, id_row, id_column);
		idrecv = id_from_index(id_row + 1, id_column); //if (idrecv == -1) {idrecv = id_from_index(id_row, 0);
		MPI_Irecv(bottom_edge_recv, N, MPI_BYTE, idrecv, tag_num = 3, MPI_COMM_WORLD, &request[2]);

		//4. send bottom
		id_to_index(id, id_row, id_column);
		idsend = id_from_index(id_row + 1, id_column); // if (idsend == -1) {idsend = id_from_index(id_row, columns - 1);}
		if (periodic == true || test == 1) { MPI_Isend(bottom_edge, N, MPI_BYTE, idsend, tag_num = 4, MPI_COMM_WORLD, &request[3]); }
		else { MPI_Isend(zeros, N, MPI_BYTE, idsend, tag_num = 4, MPI_COMM_WORLD, &request[3]); }
		//receive from top
		id_to_index(id, id_row, id_column);
		idrecv = id_from_index(id_row - 1, id_column); //if (idrecv == -1) {idrecv = id_from_index(id_row, 0);
		MPI_Irecv(top_edge_recv, N, MPI_BYTE, idrecv, tag_num = 4, MPI_COMM_WORLD, &request[3]);

		//5. send bottom right edge
		id_to_index(id, id_row, id_column);
		idsend = id_from_index(id_row + 1, id_column + 1); // if (idsend == -1) {idsend = id_from_index(id_row, columns - 1);}
		if (periodic == true || test == 1) { MPI_Isend(right_bottom, 1, MPI_BYTE, idsend, tag_num = 5, MPI_COMM_WORLD, &request[4]); }
		else { MPI_Isend(zero, 1, MPI_BYTE, idsend, tag_num = 5, MPI_COMM_WORLD, &request[4]); }
		//receive from top left
		id_to_index(id, id_row, id_column);
		idrecv = id_from_index(id_row - 1, id_column - 1); //if (idrecv == -1) {idrecv = id_from_index(id_row, 0);
		MPI_Irecv(left_top_recv, 1, MPI_BYTE, idrecv, tag_num = 5, MPI_COMM_WORLD, &request[4]);

		//6. send top left edge
		id_to_index(id, id_row, id_column);
		idsend = id_from_index(id_row - 1, id_column - 1); // if (idsend == -1) {idsend = id_from_index(id_row, columns - 1);}
		if (periodic == true || test == 1) { MPI_Isend(left_top, 1, MPI_BYTE, idsend, tag_num = 6, MPI_COMM_WORLD, &request[5]); }
		else { MPI_Isend(zero, 1, MPI_BYTE, idsend, tag_num = 6, MPI_COMM_WORLD, &request[5]); }
		//receive from bottom right edge
		id_to_index(id, id_row, id_column);
		idrecv = id_from_index(id_row + 1, id_column + 1); //if (idrecv == -1) {idrecv = id_from_index(id_row, 0);
		MPI_Irecv(right_bottom_recv, 1, MPI_BYTE, idrecv, tag_num = 6, MPI_COMM_WORLD, &request[5]);

		//7. send top right edge
		id_to_index(id, id_row, id_column);
		idsend = id_from_index(id_row - 1, id_column + 1); // if (idsend == -1) {idsend = id_from_index(id_row, columns - 1);}
		if (periodic == true || test == 1) { MPI_Isend(right_top, 1, MPI_BYTE, idsend, tag_num = 7, MPI_COMM_WORLD, &request[6]); }
		else { MPI_Isend(zero, 1, MPI_BYTE, idsend, tag_num = 7, MPI_COMM_WORLD, &request[6]); }
		//receive from bottom left edge
		id_to_index(id, id_row, id_column);
		idrecv = id_from_index(id_row + 1, id_column - 1); //if (idrecv == -1) {idrecv = id_from_index(id_row, 0);
		MPI_Irecv(left_bottom_recv, 1, MPI_BYTE, idrecv, tag_num = 7, MPI_COMM_WORLD, &request[6]);

		//8. send bottom left edge
		id_to_index(id, id_row, id_column);
		idsend = id_from_index(id_row + 1, id_column - 1); // if (idsend == -1) {idsend = id_from_index(id_row, columns - 1);}
		if (periodic == true || test == 1) { MPI_Isend(left_bottom, 1, MPI_BYTE, idsend, tag_num = 8, MPI_COMM_WORLD, &request[7]); }
		else { MPI_Isend(zero, 1, MPI_BYTE, idsend, tag_num = 8, MPI_COMM_WORLD, &request[7]); }
		//receive from top right edge
		id_to_index(id, id_row, id_column);
		idrecv = id_from_index(id_row - 1, id_column + 1); //if (idrecv == -1) {idrecv = id_from_index(id_row, 0);
		MPI_Irecv(right_top_recv, 1, MPI_BYTE, idrecv, tag_num = 8, MPI_COMM_WORLD, &request[7]);
		//----------

		MPI_Waitall(8, request, MPI_STATUS_IGNORE);

		// Fill in the boundaries of the Grid
		//-----
		for (int i = 0; i < N; i++) { grid[0][i + 1] = top_edge_recv[i]; }
		for (int i = 0; i < N; i++) { grid[N+1][i + 1] = bottom_edge_recv[i]; }

		for (int i = 0; i < N; i++) { grid[i + 1][0] = left_edge_recv[i]; }
		for (int i = 0; i < N; i++) { grid[i + 1][N+1] = right_edge_recv[i]; }

		grid[0][0] = left_top_recv[0]; 	    grid[0][N + 1] = right_top_recv[0];
		grid[N + 1][0] = left_bottom_recv[0]; 	grid[N + 1][N + 1] = right_bottom_recv[0];
		//-----

		// do iteration...
		//----------------------------------------------------------
		//cout << "it: " << n << endl;
		do_iteration();
		id_to_index(id, id_row, id_column);
		grid_to_file(n, id_row, id_column); // write output grid to a file
		//----------------------------------------------------------

		// delete the static memory
		delete[] left_edge, right_edge, right_edge_recv, left_edge_recv, top_edge, bottom_edge, bottom_edge_recv, top_edge_recv, zeros, zero;
		delete[] left_top, left_bottom, right_top, right_bottom, left_top_recv, left_bottom_recv, right_top_recv, right_bottom_recv;

	}

	MPI_Finalize();
}
//----------------------------------------------------------------