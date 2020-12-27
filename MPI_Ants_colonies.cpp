// MPI_Ants_colonies.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//
#include <algorithm>
#include <iostream>
#include <numeric>
#include "utils.h"
#include "mpi.h"

int main(int argc, char **argv)
{
	int rank, num;
	double time;

	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num);

	if (rank == 0)
	{
		printf("NbOfNodes %d\n", num);
	}

	int i, j, ant_counter, cities_counter;

	long external_loop_counter = 0;
	int *matrixTime = NULL;
	long double *tau;
	long double *eta;
	long double *tauUpdate;

	int *bestPath;
	int *otherBestPath;
	int *currentPath;
	long bestCost = MAX;
	long otherBestCost;
	long double *otherTau;
	long double *tempTau;

	long *randomNumbers = NULL;
	long nRandomNumbers = 0;
	long random_counter = 0;

	int nAnts;
	int totalNAnts;
	long externalIterations;
	long onNodeIteration;
	long double alpha;
	long double beta;
	long double evaporationCoeff;

	int *nAntsPerNode = (int *)malloc(num * sizeof(int));

	long terminationCondition = 0;

	long otherTerminationCondition = 0;
	double terminationConditionPercentage = 0.7;

	int nTasks;
	int mMachines;
	int Q;

	if (rank == 0)
	{
		nTasks = 20;
		mMachines = 5;

		matrixTime = (int *)malloc(nTasks * mMachines * sizeof(int));

		std::ifstream fin("Matr2.txt");
		for (int i = 0; i < nTasks * mMachines; i++)
		{
			int a;
			fin >> a;
			matrixTime[i] = a;
		}

		std::ifstream in;
		in.open("random.txt");
		if (!in.is_open())
		{
			printf("Cannot open file.\n");
			printf("The filepath %s is incorrect\n", "random.txt");
			MPI_Finalize();
			return -1;
		}

		char out[20];

		// Define number of random numbers
		in >> out;
		nRandomNumbers = atol(out);

		// Allocation of random numbers vector
		randomNumbers = (long *)malloc(nRandomNumbers * sizeof(long));

		for (i = 0; i < nRandomNumbers; i++)
		{
			randomNumbers[i] = 0;
		}

		i = 0;
		while (!in.eof())
		{
			in >> out;
			randomNumbers[i] = atol(out);
			i++;
		}

		in.close();

		totalNAnts = 32;
		onNodeIteration = 4000;
		externalIterations = 2;
		alpha = 2;
		beta = 2;
		evaporationCoeff = 0.0021;
	}
	/******** START TIMER ********/
	time = MPI_Wtime();
	/*****************************/

	// Share number of random numbers
	if (MPI_Bcast(&nRandomNumbers, 1, MPI_LONG, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
	{
		printf("Node %d : Error in Broadcast of nRandomNumbers", rank);
		MPI_Finalize();
		return -1;
	}
	if (MPI_Bcast(&nTasks, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
	{
		printf("Node %d : Error in Broadcast of nTasks", rank);
		MPI_Finalize();
		return -1;
	}
	if (MPI_Bcast(&mMachines, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
	{
		printf("Node %d : Error in Broadcast of mMachines", rank);
		MPI_Finalize();
		return -1;
	}
	if (MPI_Bcast(&totalNAnts, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
	{
		printf("Node %d : Error in Broadcast of totalNAnts", rank);
		MPI_Finalize();
		return -1;
	}
	if (MPI_Bcast(&onNodeIteration, 1, MPI_LONG, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
	{
		printf("Node %d : Error in Broadcast of onNodeIteration", rank);
		MPI_Finalize();
		return -1;
	}
	if (MPI_Bcast(&externalIterations, 1, MPI_LONG, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
	{
		printf("Node %d : Error in Broadcast of externalIterations", rank);
		MPI_Finalize();
		return -1;
	}
	if (MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
	{
		printf("Node %d : Error in Broadcast of alpha", rank);
		MPI_Finalize();
		return -1;
	}
	if (MPI_Bcast(&beta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
	{
		printf("Node %d : Error in Broadcast of beta", rank);
		MPI_Finalize();
		return -1;
	}
	if (MPI_Bcast(&evaporationCoeff, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
	{
		printf("Node %d : Error in Broadcast of evaporationCoeff", rank);
		MPI_Finalize();
		return -1;
	}

	// Allocation of matrixTime for non-root nodes
	if (rank != 0)
	{
		matrixTime = (int *)malloc(nTasks * mMachines * sizeof(int));
		randomNumbers = (long *)malloc(nRandomNumbers * sizeof(long));
		for (i = 0; i < nRandomNumbers; i++)
		{
			randomNumbers[i] = 0;
		}
	}

	if (MPI_Bcast(&matrixTime[0], nTasks * mMachines, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
	{
		printf("Node %d : Error in Broadcast of matrixTime", rank);
		MPI_Finalize();
		return -1;
	}

	// Share all random numbers
	if (MPI_Bcast(&randomNumbers[0], 20, MPI_LONG, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
	{
		printf("Node %d : Error in Broadcast of randomNumbers", rank);
		MPI_Finalize();
		return -1;
	}

	eta = (long double *)malloc(nTasks * nTasks * sizeof(long double));

	tau = (long double *)malloc(nTasks * nTasks * sizeof(long double));
	tauUpdate = (long double *)malloc(nTasks * nTasks * sizeof(long double));

	otherBestPath = (int *)malloc(nTasks * sizeof(int));
	otherTau = (long double *)malloc(nTasks * nTasks * sizeof(long double));
	tempTau = (long double *)malloc(nTasks * nTasks * sizeof(long double));
	for (i = 0; i < nTasks; i++)
	{
		otherBestPath[i] = -1;
	}

	bestPath = (int *)malloc(nTasks * sizeof(int));
	currentPath = (int *)malloc(nTasks * sizeof(int));

	for (i = 0; i < nTasks; i++)
	{
		currentPath[i] = -1;
		bestPath[i] = -1;
	}

	initTau(tau, nTasks);
	initEta(eta, matrixTime, nTasks, mMachines);

	Q = getLowerBound(matrixTime, nTasks, mMachines);

	int antsPerNode = totalNAnts / num;
	int restAnts = totalNAnts - antsPerNode * num;

	for (i = 0; i < num; i++)
	{
		nAntsPerNode[i] = antsPerNode;
		if (restAnts > i)
		{
			nAntsPerNode[i]++;
		}
	}

	nAnts = nAntsPerNode[rank];

	int nAntsBeforeMe = 0;
	for (i = 0; i < num; i++)
	{
		if (i < rank)
		{
			nAntsBeforeMe += nAntsPerNode[i];
		}
		else
		{
			i = num;
		}
	}

	random_counter = (random_counter + (onNodeIteration * nAntsBeforeMe * nTasks)) % nRandomNumbers;

	int antsBestCost = MAX;

	for (long external_loop_counter = 0; external_loop_counter < externalIterations; ++external_loop_counter) //&& (terminationCondition < (long)ceilf(externalIterations * onNodeIteration)))
	{

		for (long loop_counter = 0; loop_counter < onNodeIteration; loop_counter++)
		{
			for (ant_counter = 0; ant_counter < nAnts; ant_counter++)
			{
				for (i = 0; i < nTasks; i++)
				{
					currentPath[i] = -1;
				}

				//long rand = randomNumbers[random_counter];
				int currentTask = rand() % nTasks;
				//random_counter = (random_counter + 1) % nRandomNumbers;

				currentPath[currentTask] = 0;
				for (cities_counter = 1; cities_counter < nTasks; cities_counter++)
				{
					// rand = randomNumbers[random_counter];
					// random_counter = (random_counter + 1) % nRandomNumbers;
					currentTask = computeNextTask(currentTask, currentPath, tau, eta, nTasks, alpha, beta, /* rand*/ 0);

					if (currentTask == -1)
					{
						printf("There is an error choosing the next task in iteration %ld for ant %d on node %d\n", loop_counter, ant_counter, rank);
						MPI_Finalize();
						return -1;
					}

					currentPath[currentTask] = cities_counter;
				}

				long oldCost = bestCost;
				bestCost = computeCost(bestCost, bestPath, currentPath, matrixTime, nTasks, mMachines);

				if (oldCost > bestCost)
				{
					copyVectorInt(currentPath, bestPath, nTasks);
				}
			}

			if (bestCost < antsBestCost)
			{
				antsBestCost = bestCost;
				terminationCondition = 0;
			}
			else
			{
				terminationCondition++;
			}

			for (j = 0; j < nTasks * nTasks; j++)
			{
				tau[j] *= (1 - evaporationCoeff);
			}

			updateTau(tau, bestPath, bestCost, nTasks, Q);
		}

		for (int i = 0; i < nTasks; i++)
		{
			for (int j = 0; j < nTasks; ++j)
			{
				tauUpdate[j * nTasks + i] = 1.0;
				tempTau[j * nTasks + i] = 0.0;
			}
		}

		long tempBestCost = bestCost;
		int *tempBestPath = (int *)malloc(nTasks * sizeof(int));
		long tempTerminationCondition = terminationCondition;
		copyVectorInt(bestPath, tempBestPath, nTasks);

		for (i = 1; i < num; i++)
		{
			if (rank == i)
			{
				copyVectorInt(bestPath, otherBestPath, nTasks);
				otherTerminationCondition = terminationCondition;
				otherBestCost = bestCost;
				copyVectorLongdouble(tau, otherTau, nTasks * nTasks);
			}

			if (MPI_Bcast(&otherBestPath[0], nTasks, MPI_INT, i, MPI_COMM_WORLD) != MPI_SUCCESS)
			{
				printf("Node %d : Error in Broadcast of otherBestPath", rank);
				MPI_Finalize();
				return -1;
			}
			if (MPI_Bcast(&otherTau[0], nTasks * nTasks, MPI_LONG_DOUBLE, i, MPI_COMM_WORLD) != MPI_SUCCESS)
			{
				printf("Node %d : Error in Broadcast of otherTau", rank);
				MPI_Finalize();
				return -1;
			}
			if (MPI_Bcast(&otherTerminationCondition, 1, MPI_LONG, i, MPI_COMM_WORLD) != MPI_SUCCESS)
			{
				printf("Node %d : Error in Broadcast of otherTerminationCondition", rank);
				MPI_Finalize();
				return -1;
			}
			if (MPI_Bcast(&otherBestCost, 1, MPI_LONG, i, MPI_COMM_WORLD) != MPI_SUCCESS)
			{
				printf("Node %d : Error in Broadcast of otherBestCost", rank);
				MPI_Finalize();
				return -1;
			}

			// If i am not node i, I will check if values from node i are better than mine
			if (rank != i)
			{
				if (otherBestCost < tempBestCost)
				{
					tempTerminationCondition = otherTerminationCondition;
					tempBestCost = otherBestCost;
					copyVectorInt(otherBestPath, tempBestPath, nTasks);
				}
				else if (otherBestCost == tempBestCost)
				{
					tempTerminationCondition += otherTerminationCondition;
				}

				// Update pheromons received from other node
				/*for (int i = 0; i < nTasks; i++)
		 		 {
		 		 	for (int j = 0; j < nTasks; ++j)
		 		 	{
		 	 			tauUpdate[j * nTasks + i] += 1;
		 		 		tempTau[j * nTasks + i] += otherTau[j * nTasks + i];
		 		 	}
		 		 }*/
			}
		}

		// Compute the average for each pheromons value received
		for (int i = 0; i < nTasks; i++)
		{
			for (int j = 0; j < nTasks; ++j)
			{
				tau[j * nTasks + i] = /*tempTau[j * nTasks + i];*/ std::max(otherTau[j * nTasks + i], tau[j * nTasks + i]);
				// tau[j * nTasks + i] = tau[j * nTasks + i] / tauUpdate[j * nTasks + i];
			}
		}

		bestCost = tempBestCost;
		copyVectorInt(tempBestPath, bestPath, nTasks);
		terminationCondition = tempTerminationCondition;
	}

	// Merge solution into root
	if (rank == 0)
	{
		for (i = 1; i < num; i++)
		{
			if (MPI_Recv(&otherBestPath[0], nTasks, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status) != MPI_SUCCESS)
			{
				printf("Node %d : Error in Recv of otherBestPath", rank);
				MPI_Finalize();
				return -1;
			}
			long oldCost = bestCost;
			bestCost = computeCost(bestCost, bestPath, otherBestPath, matrixTime, nTasks, mMachines);

			if (oldCost > bestCost)
			{
				copyVectorInt(otherBestPath, bestPath, nTasks);
				updateTau(tau, bestPath, bestCost, nTasks, Q);
			}
		}
	}
	else
	{
		if (MPI_Send(&bestPath[0], nTasks, MPI_INT, 0, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
		{
			printf("Node %d : Error in Send of bestPath", rank);
			MPI_Finalize();
			return -1;
		}
	}

	if (rank == 0)
	{
		std::cout << "\nbestCost----->";
		std::cout << bestCost << "\n\n";
		for (int i = 0; i < nTasks; ++i)
		{
			std::cout << bestPath[i] << " ";
		}
		printf("\nTime: %f sec\n", MPI_Wtime() - time);
	}
	MPI_Finalize();
}
