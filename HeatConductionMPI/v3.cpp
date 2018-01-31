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
	double t1, t2, t3, t4, t5, t6, t9, t10, t11;
	double t7a, t7b, t7ba, t7bb, t7babis, t7za, t7zb, t7zc, t7da, t7db, t7dc, t7dd;
	double t8a, t8ba, t8bb, t8bc, t8babis, t8za, t8zb, t8fa, t8fb, t8fc, t8fd;

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
		problem.u_nplus1 = std::vector<double>(problem.s + 2);
		problem.u_n = std::vector<double>(problem.s + 2);

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
	MPI_Bcast(&problem.Text_0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	t3 = MPI_Wtime();
	tcomm += (t3-t2);
	
	// we have s and can initialize vector for the problems
	if (rank != 0){
		problem.u_nplus1 = std::vector<double>(problem.s + 2);
		problem.u_n = std::vector<double>(problem.s + 2);
	}
	
	// define final array which will store all the values and then gather
	double *final_u_n;
	final_u_n = (double*) malloc(sizeof(double) * (problem.s + 1));
	
	t4 = MPI_Wtime();
	tcomput += (t4-t3);
	
	// broadcast n = 0;
	MPI_Bcast(&problem.u_n[0], (problem.s + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	t5 = MPI_Wtime();
	tcomm += (t5-t4);
	
	// initialize linear system
	double m = 0;
	std::vector<double> a = std::vector<double>(problem.s + 1);
	std::vector<double> b = std::vector<double>(problem.s + 1);
	std::vector<double> c = std::vector<double>(problem.s + 1);
	std::vector<double> d = std::vector<double>(problem.s + 1);
	for (int i = 0; i < problem.s - 1; i++){
		a[i] = -problem.r/2;
		b[i] = problem.r + 1;
		c[i] = -problem.r/2;
		d[i] = (problem.r/2)*problem.u_n[i+2] + (1-problem.r)*problem.u_n[i + 1] + (problem.r/2)*problem.u_n[i];
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
	
	t6 = MPI_Wtime();
	tcomput += (t6-t5);
	
	/*************** Calculte n = 1 and so on ***************/
	for (int j = 1; j < (problem.n + 1); j++){
		t7a = MPI_Wtime();
		
		//Boundaries conditions
		d[0] += problem.Text_0 * (problem.r/2);
		d[problem.s - 2] += problem.Text_0 * (problem.r/2);
				
		// for a n, we compute n+1
		//forward phase
		if (rank == 0) {	
			if (npes > 1){
				for (int k = beginning + 1; k < ending - 1; k++){
					m = a[k] / b[k - 1];
					b[k] = b[k] - (m*c[k - 1]);
					d[k] = d[k] - (m*d[k - 1]);
				}
				
				t7ba = MPI_Wtime();
				tcomput += (t7ba-t7a);
				
				MPI_Send(&b[ending - 2], 1, MPI_DOUBLE, rank + 1, 100, MPI_COMM_WORLD);
				MPI_Send(&d[ending - 2], 1, MPI_DOUBLE, rank + 1, 102, MPI_COMM_WORLD);
				
				t7bb = MPI_Wtime();
				tcomm += (t7bb-t7ba);
			}
			else {
				for (int k = beginning + 1; k < ending - 2; k++){
					m = a[k] / b[k - 1];
					b[k] = b[k] - (m*c[k - 1]);
					d[k] = d[k] - (m*d[k - 1]);
				}
				
				t7babis = MPI_Wtime();
				tcomput += (t7babis-t7a);
			}
		}
		else if (rank == (npes - 1) && (rank != 0)) {
			t7za = MPI_Wtime();
			tcomput += (t7za-t7a);
				
			MPI_Recv(&b[beginning - 2], 1, MPI_DOUBLE, rank - 1, 100, MPI_COMM_WORLD, &status);
			MPI_Recv(&d[beginning - 2], 1, MPI_DOUBLE, rank - 1, 102, MPI_COMM_WORLD, &status);
			
			t7zb = MPI_Wtime();
			tcomm += (t7zb-t7za);
				
			for (int k = beginning - 1; k < ending - 2; k++){
				m = a[k] / b[k - 1];
				b[k] = b[k] - (m*c[k - 1]);
				d[k] = d[k] - (m*d[k - 1]);
			}
			
			t7zc = MPI_Wtime();
			tcomput += (t7zc-t7zb);
		}
		else {
			t7da = MPI_Wtime();
			tcomput += (t7da-t7a);
			
			MPI_Recv(&b[beginning - 2], 1, MPI_DOUBLE, rank - 1, 100, MPI_COMM_WORLD, &status);
			MPI_Recv(&d[beginning - 2], 1, MPI_DOUBLE, rank - 1, 102, MPI_COMM_WORLD, &status);
			
			t7db = MPI_Wtime();
			tcomm += (t7db-t7da);
					
			for (int k = beginning - 1; k < ending - 1; k++){
				m = a[k] / b[k - 1];
				b[k] = b[k] - (m*c[k - 1]);
				d[k] = d[k] - (m*d[k - 1]);
			}
			
			t7dc = MPI_Wtime();
			tcomput += (t7dc-t7db);
			
			MPI_Send(&b[ending - 2], 1, MPI_DOUBLE, rank + 1, 100, MPI_COMM_WORLD);
			MPI_Send(&d[ending - 2], 1, MPI_DOUBLE, rank + 1, 102, MPI_COMM_WORLD);
			
			t7dd = MPI_Wtime();
			tcomm += (t7dd-t7dc);
		}

		t8a = MPI_Wtime();

		//backward phase
		problem.u_nplus1[problem.s] = problem.Text_0;
		problem.u_nplus1[0] = problem.Text_0;
		problem.u_nplus1[problem.s - 1] = d[problem.s - 2] / b[problem.s - 2]; 
		if (rank == 0) {	
			if (npes > 1){
				t8ba = MPI_Wtime();
				tcomput += (t8ba-t8a);
			
				MPI_Recv(&problem.u_nplus1[ending], 1, MPI_DOUBLE, rank + 1, 200, MPI_COMM_WORLD, &status);
				
				t8bb = MPI_Wtime();
				tcomm += (t8bb-t8ba);
				
				for (int k = ending - 2; k > beginning - 1; k--){
					problem.u_nplus1[k + 1] = (d[k] - (c[k] * problem.u_nplus1[k + 2])) / b[k];
				}
				
				t8bc = MPI_Wtime();
				tcomput += (t8bc-t8bb);
			}
			else {
				for (int k = ending - 4; k > beginning - 1; k--){
					problem.u_nplus1[k + 1] = (d[k] - (c[k] * problem.u_nplus1[k + 2])) / b[k];
				}
				
				t8babis = MPI_Wtime();
				tcomput += (t8babis-t8a);
			}			
		}
		else if (rank == (npes - 1) && (rank != 0)) {
			for (int k = ending - 4; k > beginning - 2; k--){
				problem.u_nplus1[k + 1] = (d[k] - (c[k] * problem.u_nplus1[k + 2])) / b[k];
			}
			
			t8za = MPI_Wtime();
			tcomput += (t8za-t8a);
			
			MPI_Send(&problem.u_nplus1[beginning], 1, MPI_DOUBLE, rank - 1, 200, MPI_COMM_WORLD);
			
			t8zb = MPI_Wtime();
			tcomm += (t8zb-t8za);	
		}
		else {
			t8fa = MPI_Wtime();
			tcomput += (t8fa-t8a);
			
			MPI_Recv(&problem.u_nplus1[ending], 1, MPI_DOUBLE, rank + 1, 200, MPI_COMM_WORLD, &status);
			
			t8fb = MPI_Wtime();
			tcomm += (t8fb-t8fa);
						
			for (int k = ending - 2; k > beginning - 2; k--){
				problem.u_nplus1[k + 1] = (d[k] - (c[k] * problem.u_nplus1[k + 2])) / b[k];
			}
			
			t8fc = MPI_Wtime();
			tcomput += (t8fc-t8fb);
			
			MPI_Send(&problem.u_nplus1[beginning], 1, MPI_DOUBLE, rank - 1, 200, MPI_COMM_WORLD);
			
			t8fd = MPI_Wtime();
			tcomm += (t8fd-t8fc);	
		}

		t9 = MPI_Wtime();

		// we can update u_n
		problem.u_n = problem.u_nplus1;

		t10 = MPI_Wtime();
		tcomput += (t10-t9);

		/* Process 0 gather data from every other process */
		// need to relocate computed data at the beginning of each u_n
		std::vector<double> temp(problem.u_n.begin() + beginning, problem.u_n.begin() + ending);
		std::copy(temp.begin(), temp.end(), problem.u_n.begin());
	
		// gather the data
		MPI_Gatherv(&problem.u_n[0], sendcounts[rank], MPI_DOUBLE, final_u_n, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		// broadcast n once back into a vector;
		if (rank == 0) {
			for (int i = 0; i < problem.s + 1; i++) problem.u_n[i] = final_u_n[i];
		}
		MPI_Bcast(&problem.u_n[0], (problem.s + 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
		
		t10 = MPI_Wtime();
		tcomm += (t10-t9);
		
		// we set back correctly vector b and d
		for (int i = 0; i < problem.s - 1; i++){
			b[i] = problem.r + 1; // central diagonal
			d[i] = (problem.r/2)*problem.u_n[i+2] + (1-problem.r)*problem.u_n[i + 1] + (problem.r/2)*problem.u_n[i];
		}
		
		t11 = MPI_Wtime();
		tcomput += (t11-t10);
	}

	/*************** Process 0 write properly the results in a file ***************/
	if (rank == 0) {
		FILE *fileFTCS; // file to store data
		std::string name = "CRANK-";
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

	printf("Elapsed time: %f, rank: %d\n", t11-t1, rank);
	printf("Tcomm: %f  rank: %d\n", tcomm, rank);
	printf("Tcomput: %f  rank: %d\n", tcomput, rank);
	
	/*************** Terminate the program ***************/
	free(sendcounts);
	free(displs);
	free(final_u_n);
	MPI_Finalize();
	return 0;
}
