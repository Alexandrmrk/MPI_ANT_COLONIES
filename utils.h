#pragma once
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <limits>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>

const long double INF = std::numeric_limits<long double>::infinity();
const int MAX = std::numeric_limits<int>::max();

int getMatrixIndex(int i, int j, int matrixSize)
{
	return (i * matrixSize) + j;
}

void copyVectorInt(int *in, int *out, int size)
{
	int k;
	for (k = 0; k < size; k++)
	{
		out[k] = (int)in[k];
	}
}

void copyVectordouble(double *in, double *out, int size)
{
	int k;
	for (k = 0; k < size; k++)
	{
		out[k] = (double)in[k];
	}
}

void copyVectorLongdouble(long double *in, long double *out, int size)
{
	int k;
	for (k = 0; k < size; k++)
	{
		out[k] = (long double)in[k];
	}
}

void computeProbabilities(int currentTask, long double *probabilities, int *path, long double *tau, long double *eta, int nTasks, long double alpha, long double beta)
{
	int i;
	long double total = 0;
	for (i = 0; i < nTasks; i++)
	{
		if (path[i] != -1 || i == currentTask)
		{
			probabilities[i] = 0.0;
		}
		else
		{
			long double p = pow(tau[getMatrixIndex(currentTask, i, nTasks)], alpha) * pow(eta[getMatrixIndex(currentTask, i, nTasks)], beta);
			probabilities[i] = p;
			total += p;
		}
	}

	if (total == 0)
	{
		i = 0;
		for (i = 0; i < nTasks; i++)
		{
			if (path[i] == -1 && i != currentTask)
			{
				probabilities[i] = 1.0;
				total++;
			}
		}
	}

	for (i = 0; i < nTasks; i++)
	{
		probabilities[i] = probabilities[i] / total;
	}
}
int computeNextTask(int currentTask, int *path, long double *tau, long double *eta, int nTasks, long double alpha, long double beta, long random)
{
	int i = 0;
	long double *probabilities;
	probabilities = (long double *)malloc(nTasks * sizeof(long double));
	computeProbabilities(currentTask, probabilities, path, tau, eta, nTasks, alpha, beta);

	int value = (rand() % 100) + 1;
	long double sum = 0;

	for (i = 0; i < nTasks; i++)
	{
		sum += ceilf(probabilities[i] * 100);
		if (sum >= value)
		{
			free(probabilities);
			return i;
		}
	}
	free(probabilities);
	return -1;
}

long computeCost(long bestCost, int *bestPath, int *currentPath, int *matrixTime, int nTasks, int mMachines)
{
	int i;
	long currentCost = 0;
	int *matrixX;
	matrixX = (int *)malloc(nTasks * mMachines * sizeof(int));

	int *completedTasks = (int *)malloc(nTasks * sizeof(int));

	for (i = 0; i < nTasks; i++)
	{
		completedTasks[currentPath[i]] = i;
	}

	for (int i = 0; i < nTasks * mMachines; ++i)
	{
		matrixX[i] = 0;
	}

	for (int i = 1; i < nTasks; ++i)
	{
		int val = matrixX[i - 1] + matrixTime[completedTasks[i - 1]];
		matrixX[i] = val;
	}

	for (int j = 1; j < mMachines; ++j)
	{
		matrixX[j * nTasks] = matrixX[(j - 1) * nTasks + 1];
		for (int i = 1; i < nTasks; ++i)
		{
			int a = matrixX[(j - 1) * nTasks + i] + matrixTime[(j - 1) * nTasks + completedTasks[i]];
			int b = matrixX[j * nTasks + i - 1] + matrixTime[j * nTasks + completedTasks[i - 1]];
			if (a > b)
			{
				matrixX[j * nTasks + i] = a;
			}
			else
			{
				matrixX[j * nTasks + i] = b;
			}
		}
	}

	currentCost = matrixX[(mMachines - 1) * nTasks + nTasks - 1] + matrixTime[(mMachines - 1) * nTasks + completedTasks[nTasks - 1]];

	free(matrixX);
	free(completedTasks);

	if (bestCost > currentCost)
	{
		return currentCost;
	}
	else
	{
		return bestCost;
	}
}

void updateTau(long double *tau, int *path, long cost, int nTasks, int Q)
{
	int i;
	int *completedTasks = (int *)malloc(nTasks * sizeof(int));

	for (i = 0; i < nTasks; i++)
	{
		completedTasks[path[i]] = i;
	}

	for (int i = 0; i < nTasks - 1; ++i)
	{
		tau[getMatrixIndex(completedTasks[i], completedTasks[i + 1], nTasks)] += (long double)Q / cost;
	}
	for (int i = 0; i < nTasks - 1; ++i)
	{
		if (tau[getMatrixIndex(completedTasks[i], completedTasks[i + 1], nTasks)] > 1)
		{
			tau[getMatrixIndex(completedTasks[i], completedTasks[i + 1], nTasks)] = 1.0;
		}
	}
	free(completedTasks);
}

void initTau(long double *tau, int nTask)
{
	for (int i = 0; i < nTask; ++i)
	{
		for (int j = 0; j < nTask; ++j)
		{
			tau[i * nTask + j] = 1.0 / nTask;
		}
	}
}

void initEta(long double *eta, int *matrixTime, int nTasks, int mMachines)
{
	std::vector<std::vector<int>> mT;
	mT.resize(nTasks);
	for (int i = 0; i < nTasks; ++i)
		mT[i].resize(mMachines);

	for (int i = 0; i < nTasks; i++)
	{
		for (int j = 0; j < mMachines; ++j)
		{
			mT[i][j] = matrixTime[j * nTasks + i];
		}
	}

	for (int i = 0; i < nTasks; ++i)
	{
		for (int j = 0; j < nTasks; ++j)
		{
			if (i != j)
			{
				int sum = 0;
				int a = mT[i][0];
				int b = a;
				for (int k = 1; k < mMachines; ++k)
				{
					a += mT[i][k];
					b += mT[j][k - 1];
					sum += std::abs(a - b);
					if (a > b)
					{
						b = a;
					}
				}
				eta[getMatrixIndex(i, j, nTasks)] = (long double)1 / sum;
			}
			else
			{
				eta[getMatrixIndex(i, j, nTasks)] = INF;
			}
		}
	}
}

int getLowerBound(int *matrixTime, int nTask, int mMachines)
{
	int sum = 0;
	int min = MAX;
	for (int i = 0; i < nTask; ++i)
	{
		for (int j = 0; j < mMachines; ++j)
		{
			if (min > matrixTime[j * nTask + i])
			{
				min = matrixTime[j * nTask + i];
			}
		}
	}
	sum = min * (mMachines - 1);

	for (int i = 0; i < nTask; ++i)
	{
		sum += matrixTime[(mMachines - 1) * nTask + i];
	}
	return sum;
}
