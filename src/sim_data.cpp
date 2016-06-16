/**
 * 
	The MIT License (MIT)

	Copyright (c) 2016 Dylan

	Permission is hereby granted, free of charge, to any person obtaining a copy
	of this software and associated documentation files (the "Software"), to deal
	in the Software without restriction, including without limitation the rights
	to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
	copies of the Software, and to permit persons to whom the Software is
	furnished to do so, subject to the following conditions:

	The above copyright notice and this permission notice shall be included in all
	copies or substantial portions of the Software.

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
	OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
	SOFTWARE.
*/

#include "../include/sim_data.h"

#include<stdlib.h>
#include<iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <omp.h>

#include "mkl.h"


SimulationData::SimulationData(int N, int M, double U, double J) {

	this->N = N;
	this->M = M;
	this->U = U;
	this->J = J;
	this->D = floor(tgamma(N + M) / (tgamma(N+1)*tgamma(M)));

	this->eigenvalues = (double*)mkl_malloc(D*sizeof(double), 64);
	this->eigenvectors = (double*)mkl_malloc(D*sizeof(double), 64);
	this->ground_state = (double*)mkl_malloc(D*sizeof(double), 64);

	









}

SimulationData::~SimulationData(){

	mkl_free(eigenvectors);
	mkl_free(eigenvalues);
	mkl_free(ground_state);
}

