//#include "C:\Program Files (x86)\Microsoft SDKs\MPI\Include\mpi.h"
//#include "mpi.h"

#include <vector>
#include <stdlib.h> // for malloc
#include <algorithm> // for std::copy

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

int main2(){

	int rank, npes;
	MPI_Status status;

	MPI_Init(0, 0);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &npes);

	/* each process create a problem */
	struct Problem problem;
	problem.u_nplus1 = std::vector<double>(problem.s + 1);
	problem.u_n = std::vector<double>(problem.s + 1);
	problem.sub_u_n = std::vector<double>(problem.s + 1);

	/* Process 0 get the full problem */
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
		problem.u_nplus1 = std::vector<double>(problem.s + 1);
		problem.u_n = std::vector<double>(problem.s + 1);

		/* initialisation n = 0 */
		for (int i = 1; i < problem.s; i++){
			problem.u_n[i] = Tin_0;
		}
		problem.u_n[0] = problem.Text_0;
		problem.u_n[problem.s] = problem.Text_0;
		printf("u_n %d %d\n", rank, problem.u_n[0]);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	/* Process 0 broadcast some values needed and scatter vectors */
	MPI_Bcast(&problem.n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&problem.s, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&problem.r, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&problem.u_n[0], (problem.s + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
	printf("s %d %d\n", rank, problem.s);
	printf("u_n %d %d\n", rank, problem.u_n[0]);
	int *sendcounts;    // array describing how many elements to send to each process
	int *displs; // array describing the displacements where each segment begins
	int rem = (problem.s + 1) % npes; // elements remaining after division among processes
	int sum = 0; // Sum of counts. Used to calculate displacements
	sendcounts = (int*)malloc(sizeof(int)*npes);
	displs = (int*)malloc(sizeof(int)*npes);
	// calculate send counts and displacements
	for (int i = 0; i < npes; i++) {
		sendcounts[i] = (problem.s + 1) / npes;
		if (rem > 0) {
			sendcounts[i]++;
			rem--;
		}
		displs[i] = sum;
		sum += sendcounts[i];
	}

	MPI_Scatterv(&problem.u_n[0], sendcounts, displs, MPI_DOUBLE, &problem.sub_u_n[0], (problem.s + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
	// define beginning and ending of each portions for each threads
	int beginning = displs[rank];
	int ending;
	if (rank == (npes - 1)){
		ending = sum;
	}
	else {
		ending = displs[rank + 1];
	}
	printf("%d %d %d\n", rank, beginning, ending);
	printf("sub %d %d\n", rank, problem.sub_u_n[0]);

	/* Reorganize the vector u_n */
	std::copy((problem.sub_u_n).begin(), (problem.sub_u_n).end(), ((problem.u_n).begin() + beginning));
	printf("ok %d\n", rank);
	/* Terminate the program */
	free(sendcounts);
	free(displs);
	MPI_Finalize();
	return 0;
}
