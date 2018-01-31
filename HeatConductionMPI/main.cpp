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
	double t1, t2, t3, t4, t5, t6, t7a, t7b, t7c, t7d, t8, t9, t10;

	MPI_Init(0, 0);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &npes);

	double tcomm = 0; 
	double tcomput = 0;

	t1 = MPI_Wtime();

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
		double dx = 0.005; // space step = 0.05;
		double dt = 0.001; // time step = 0.01;

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
		problem.u_nplus1 = std::vector<double>(problem.s + 1);
		problem.u_n = std::vector<double>(problem.s + 1);

		// initialisation n = 0
		for (int i = 1; i < problem.s; i++){
			problem.u_n[i] = problem.Tin_0;
		}
		problem.u_n[0] = problem.Text_0;
		problem.u_n[problem.s] = problem.Text_0;
	}
	
	t2 = MPI_Wtime();
	tcomput += (t2-t1);

	/*************** Process 0 broadcast ***************/
	// broadcast values
	MPI_Bcast(&problem.n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&problem.s, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&problem.r, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	t3 = MPI_Wtime();
	tcomm += (t3-t2);
	
	// we have s and can initialize vector for the problems
	if (rank != 0){
		problem.u_nplus1 = std::vector<double>(problem.s + 1);
		problem.u_n = std::vector<double>(problem.s + 1);
	}

	t4 = MPI_Wtime();
	tcomput += (t4-t3);
	
	// broadcast n = 0;
	MPI_Bcast(&problem.u_n[0], (problem.s + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);

	t5 = MPI_Wtime();
	tcomm += (t5-t4);	

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

	t6 = MPI_Wtime();
	tcomput += (t6-t5);

	/*************** Calculte n = 1 and so on ***************/
	for (int j = 1; j < (problem.n + 1); j++){
		t7a = MPI_Wtime();
		
		// for a n, we compute n+1
		if (beginning == 0) {
			problem.u_nplus1[0] = problem.u_n[0];
			for (int i = 1; i < ending; i++){
				problem.u_nplus1[i] = problem.u_n[i] + problem.r * (problem.u_n[i + 1] - (2 * problem.u_n[i]) + problem.u_n[i - 1]); // u_nplus1 is define accrding the scheme used
			}
		}
		else {
			for (int i = beginning; i < ending; i++){
				problem.u_nplus1[i] = problem.u_n[i] + problem.r * (problem.u_n[i + 1] - (2 * problem.u_n[i]) + problem.u_n[i - 1]); // u_nplus1 is define accrding the scheme used
			}
		}
		problem.u_nplus1[problem.s] = problem.u_n[problem.s];

		t7b = MPI_Wtime();
		tcomput += (t7b-t7a);

		// we need to exchange our information with the neighbour (if there are neighboors)
		if (npes > 1) {
			// send forward, receive backward
			if ((rank != (npes - 1)) && (rank != 0)) {
				MPI_Sendrecv(&problem.u_nplus1[ending - 1], 1, MPI_DOUBLE, rank + 1, 1, &problem.u_nplus1[beginning - 1], 1, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &status);
			}
			else {
				if (rank == 0) MPI_Send(&problem.u_nplus1[ending - 1], 1, MPI_DOUBLE, (rank + 1), 1, MPI_COMM_WORLD);
				if (rank == (npes - 1)) MPI_Recv(&problem.u_nplus1[beginning - 1], 1, MPI_DOUBLE, (npes - 2), 1, MPI_COMM_WORLD, &status);
			}	
			// sent backward, receive forward
			if ((rank != (npes - 1)) && (rank != 0)) {
				MPI_Sendrecv(&problem.u_nplus1[beginning], 1, MPI_DOUBLE, rank - 1, 2, &problem.u_nplus1[ending], 1, MPI_DOUBLE, rank + 1, 2, MPI_COMM_WORLD, &status);
			}
			else {
				if (rank == 0) MPI_Recv(&problem.u_nplus1[ending], 1, MPI_DOUBLE, (rank + 1), 2, MPI_COMM_WORLD, &status);
				if (rank == (npes - 1)) MPI_Send(&problem.u_nplus1[beginning], 1, MPI_DOUBLE, (npes - 2), 2, MPI_COMM_WORLD);
			}
		}
		
		t7c = MPI_Wtime();
		tcomm += (t7c-t7b);
		
		// we can update u_n
		problem.u_n = problem.u_nplus1;
		
		t7d = MPI_Wtime();
		tcomput += (t7d-t7c);
		
	}
	
	/*************** Process 0 gather data from every other process ***************/
	
	t8 = MPI_Wtime();

	// need to relocate computed data at the beginning of each u_n	
	std::vector<int> temp(problem.u_n.begin() + beginning, problem.u_n.begin() + ending);
	std::copy(temp.begin(), temp.end(), problem.u_n.begin());
	
	t9 = MPI_Wtime();
	tcomput += (t9-t8);
	
	// define final array which will store all the values and then gather
	double *final_u_n;
	final_u_n = (double*) malloc(sizeof(double) * (problem.s + 1));
	MPI_Gatherv(&problem.u_n[0], sendcounts[rank], MPI_DOUBLE, final_u_n, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	t10 = MPI_Wtime();
	tcomm += (t10-t9);

	/*************** Process 0 write properly the results in a file ***************/
	if (rank == 0) {
		FILE *fileFTCS; // file to store data
		std::string name = "FTCS-";
		std::string dxS = std::to_string(problem.dx);
		name.append(dxS);
		char* dxC = (char*)name.c_str();
		fileFTCS = fopen( dxC , "w");
		double x = 0;
		for (int i = 0; i < (problem.s + 1); i++){
			fprintf(fileFTCS, "%f %f\n", x, final_u_n[i]);
			x += problem.dx;
		}
		fclose(fileFTCS);
	} 

	printf("Elapsed time: %f  rank: %d\n", t10-t1, rank);
	printf("Tcomm: %f  rank: %d\n", tcomm, rank);
	printf("Tcomput: %f  rank: %d\n", tcomput, rank);
	
	/*************** Terminate the program ***************/
	free(sendcounts);
	free(displs);
	free(final_u_n);
	MPI_Finalize();
	return 0;
}
