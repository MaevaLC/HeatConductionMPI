//#include "C:\Program Files (x86)\Microsoft SDKs\MPI\Include\mpi.h"
#include "mpi.h"

#include <vector>
#include <stdlib.h> // for malloc
#include <algorithm> // for std::copy
#include <string> // for std::to_string

struct Problem {
	double Tin_0; //!< initial condition Temperature
	double Text_0; //!< initial condition Temperature
	double Xmin; //!< initial condition Position
	double Xmax; //!< initial condition Position
	double Tend; //!< initial condition Time
	double D; //!< initial condition D
	double dx; //!< space step
	double dt; //!< time step
	int n; //!< number of time steps
	int s; //!< number of space steps
	double r; //!< calculation made once instead of multiple time
	std::vector<double> u_nplus1; //!< solution values vector n+1
	std::vector<double> u_n; //!< solution values vector n
	std::vector<double> sub_u_n;
};

int main(){
	int rank, npes;
	MPI_Status status;
	double t1, t2;

	MPI_Init(0, 0);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &npes);

	/*************** each process create a problem ***************/
	struct Problem problem;

	/*************** Process 0 get the full problem ***************/
	if (rank == 0) {
		double Tin_0 = 100; // initial temperature inside : 100°F
		double Text_0 = 300; // initial temperature outside : 300°F
		double Xmin = 0; // position at the left : 0 ft
		double Xmax = 1; // position at the right : 1 ft
		double Tend = 0.5; // end time of simulation : 0.5h
		double D = 0.1; // coefficient D = 0.1 ft²/h
		double dx = 0.05; // space step = 0.05;
		double dt = 0.01; // time step = 0.01;

		problem.Tin_0 = Tin_0;
		problem.Text_0 = Text_0;
		problem.Xmin = Xmin;
		problem.Xmax = Xmax;
		problem.Tend = Tend;
		problem.D = D;
		problem.dx = dx;
		problem.dt = dt;
		problem.n = int(Tend / dt);
		problem.s = int((Xmax - Xmin) / dx);
		problem.r = (D*dt) / (dx*dx);
		problem.u_nplus1 = std::vector<double>(problem.s + 2);
		problem.u_n = std::vector<double>(problem.s + 2);

		// initialisation n = 0
		for (int i = 1; i < problem.s; i++){
			problem.u_n[i] = problem.Tin_0;
		}
		problem.u_n[0] = problem.Text_0;
		problem.u_n[problem.s] = problem.Text_0;
	}

	t1 = MPI_Wtime();

	/*************** Process 0 broadcast ***************/
	// broadcast values
	MPI_Bcast(&problem.n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&problem.s, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&problem.r, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&problem.Text_0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// we have s and can initialize vector for the problems
	if (rank != 0){
		problem.u_nplus1 = std::vector<double>(problem.s + 2);
		problem.u_n = std::vector<double>(problem.s + 2);
	}

	// broadcast n = 0;
	MPI_Bcast(&problem.u_n[0], (problem.s + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// initialize linear system
	double m = 0;
	std::vector<double> a = std::vector<double>(problem.s + 1);
	std::vector<double> b = std::vector<double>(problem.s + 1);
	std::vector<double> c = std::vector<double>(problem.s + 1);
	std::vector<double> d = std::vector<double>(problem.s + 1);
	for (int i = 0; i < problem.s - 1; i++){
		a[i] = -problem.r;
		b[i] = 2 * problem.r + 1;
		c[i] = -problem.r;
		d[i] = problem.u_n[i + 1];
	}

	// calculate workload per process
	int *sendcounts;    // array describing how many elements to send to each process
	int *displs; // array describing the displacements where each segment begins
	int rem = (problem.s + 1) % npes; // elements remaining after division among processes
	int sum = 0; // Sum of counts. Used to calculate displacements
	sendcounts = (int*)malloc(sizeof(int)*npes);
	displs = (int*)malloc(sizeof(int)*npes);
	for (int i = 0; i < npes; i++) {
		sendcounts[i] = (problem.s + 1) / npes;
		if (rem > 0) {
			sendcounts[i]++;
			rem--;
		}
		displs[i] = sum;
		sum += sendcounts[i];
	}

	// define beginning and ending of each portions for each threads
	int beginning = displs[rank];
	int ending;
	if (rank == (npes - 1)){
		ending = sum;
	}
	else {
		ending = displs[rank + 1];
	}

	/*************** Calculte n = 1 and so on ***************/
	for (int j = 1; j < (problem.n + 1); j++){
		//Boundaries conditions
		d[0] += problem.Text_0 * problem.r;
		d[problem.s - 2] += problem.Text_0 * problem.r;

		// for a n, we compute n+1
		//forward phase
		if (rank == 0) {
			if (npes > 1){
				for (int k = beginning + 1; k < ending - 1; k++){
					m = a[k] / b[k - 1];
					b[k] = b[k] - (m*c[k - 1]);
					d[k] = d[k] - (m*d[k - 1]);
				}
				MPI_Send(&b[ending - 2], 1, MPI_DOUBLE, rank + 1, 100, MPI_COMM_WORLD);
				MPI_Send(&d[ending - 2], 1, MPI_DOUBLE, rank + 1, 102, MPI_COMM_WORLD);
			}
			else {
				for (int k = beginning + 1; k < ending - 2; k++){
					m = a[k] / b[k - 1];
					b[k] = b[k] - (m*c[k - 1]);
					d[k] = d[k] - (m*d[k - 1]);
				}
			}
		}
		else if (rank == (npes - 1) && (rank != 0)) {
			MPI_Recv(&b[beginning - 2], 1, MPI_DOUBLE, rank - 1, 100, MPI_COMM_WORLD, &status);
			MPI_Recv(&d[beginning - 2], 1, MPI_DOUBLE, rank - 1, 102, MPI_COMM_WORLD, &status);
			for (int k = beginning - 1; k < ending - 2; k++){
				m = a[k] / b[k - 1];
				b[k] = b[k] - (m*c[k - 1]);
				d[k] = d[k] - (m*d[k - 1]);
			}
		}
		else {
			MPI_Recv(&b[beginning - 2], 1, MPI_DOUBLE, rank - 1, 100, MPI_COMM_WORLD, &status);
			MPI_Recv(&d[beginning - 2], 1, MPI_DOUBLE, rank - 1, 102, MPI_COMM_WORLD, &status);
			for (int k = beginning - 1; k < ending - 1; k++){
				m = a[k] / b[k - 1];
				b[k] = b[k] - (m*c[k - 1]);
				d[k] = d[k] - (m*d[k - 1]);
			}
			MPI_Send(&b[ending - 2], 1, MPI_DOUBLE, rank + 1, 100, MPI_COMM_WORLD);
			MPI_Send(&d[ending - 2], 1, MPI_DOUBLE, rank + 1, 102, MPI_COMM_WORLD);
		}

		//backward phase
		problem.u_nplus1[problem.s] = problem.Text_0;
		problem.u_nplus1[0] = problem.Text_0;
		problem.u_nplus1[problem.s - 1] = d[problem.s - 2] / b[problem.s - 2];
		if (rank == 0) {
			if (npes > 1){
				MPI_Recv(&problem.u_nplus1[ending], 1, MPI_DOUBLE, rank + 1, 200, MPI_COMM_WORLD, &status);
				for (int k = ending - 2; k > beginning - 1; k--){
					problem.u_nplus1[k + 1] = (d[k] - (c[k] * problem.u_nplus1[k + 2])) / b[k];
				}
			}
			else {
				for (int k = ending - 4; k > beginning - 1; k--){
					problem.u_nplus1[k + 1] = (d[k] - (c[k] * problem.u_nplus1[k + 2])) / b[k];
				}
			}
		}
		else if (rank == (npes - 1) && (rank != 0)) {
			for (int k = ending - 4; k > beginning - 2; k--){
				problem.u_nplus1[k + 1] = (d[k] - (c[k] * problem.u_nplus1[k + 2])) / b[k];
			}
			MPI_Send(&problem.u_nplus1[beginning], 1, MPI_DOUBLE, rank - 1, 200, MPI_COMM_WORLD);
		}
		else {
			MPI_Recv(&problem.u_nplus1[ending], 1, MPI_DOUBLE, rank + 1, 200, MPI_COMM_WORLD, &status);
			for (int k = ending - 2; k > beginning - 2; k--){
				problem.u_nplus1[k + 1] = (d[k] - (c[k] * problem.u_nplus1[k + 2])) / b[k];
			}
			MPI_Send(&problem.u_nplus1[beginning], 1, MPI_DOUBLE, rank - 1, 200, MPI_COMM_WORLD);
		}

		// we can update u_n
		problem.u_n = problem.u_nplus1;

		// we set back correctly the vector b
		for (int i = 0; i < problem.s - 1; i++){
			b[i] = 2 * problem.r + 1; // central diagonal
			d[i] = problem.u_n[i + 1];
		}
	}

	/*************** Process 0 gather data from every other process ***************/
	// need to relocate computed data at the beginning of each u_n
	std::vector<int> temp(problem.u_n.begin() + beginning, problem.u_n.begin() + ending);
	std::copy(temp.begin(), temp.end(), problem.u_n.begin());

	// define final array which will store all the values and then gather
	double *final_u_n;
	final_u_n = (double*)malloc(sizeof(double)* (problem.s + 1));
	MPI_Gatherv(&problem.u_n[0], sendcounts[rank], MPI_DOUBLE, final_u_n, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	t2 = MPI_Wtime();

	/*************** Process 0 write properly the results in a file ***************/
	if (rank == 0) {
		FILE *fileFTCS; // file to store data
		std::string name = "LAAS-";
		std::string dxS = std::to_string(problem.dx);
		name.append(dxS);
		char* dxC = (char*)name.c_str();
		fileFTCS = fopen(dxC, "w");
		double x = 0;
		for (int i = 0; i < (problem.s + 1); i++){
			fprintf(fileFTCS, "%f %f\n", x, final_u_n[i]);
			x += problem.dx;
		}
		fclose(fileFTCS);
	}

	printf("Elapsed time: %f, rank: %d\n", t2 - t1, rank);

	/*************** Terminate the program ***************/
	free(sendcounts);
	free(displs);
	free(final_u_n);
	MPI_Finalize();
	return 0;
}
