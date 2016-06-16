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

#include "../include/basis.h"

#include<stdlib.h>
#include<iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <omp.h>
#include <algorithm>

#include "mkl.h"

Basis::Basis(SimulationData &simdata) {

	this->basis_vectors = (double*)mkl_malloc(simdata.D*simdata.M*sizeof(double), 64);
	this->tags = (double*)mkl_malloc(simdata.D*sizeof(double), 64);
	this->tags_sorted = (double*)mkl_malloc(simdata.D*sizeof(double), 64);


};	

void Basis::ConstructBasis(SimulationData &simdata){

	double *prev_basis;
	double *new_basis;
	prev_basis = (double*)mkl_malloc(simdata.D*simdata.N*sizeof(double), 64);
	new_basis = (double*)mkl_malloc(simdata.D*simdata.N*sizeof(double), 64);

	int index;

	this->basis_vectors[0] = simdata.N;

	//First basis vector |N,0,...,0>
	int i;
	#pragma omp parallel for private(i)
	for (i = 1; i < simdata.M; ++i) {
		this->basis_vectors[i] = 0;
	}

//	#pragma omp parallel for private(i)
//	for (i = (simdata.D - simdata.M); i < simdata.M-1; ++i) {
//		this->basis_vectors[i] = 0;
//	}
//	this->basis_vectors[simdata.D-1] = simdata.N;

	double *myPointer = 0;
	myPointer = (double*)mkl_malloc(simdata.M*sizeof(double), 64);
	double *myPointer2 = 0;
	myPointer2 = (double*)mkl_malloc(simdata.M*sizeof(double), 64);
	double *nonzero;

	int n_k;
	//Outer loop sets the basis vector
	for (int i = simdata.M; i < simdata.D*simdata.M; i+=simdata.M) {
		//Inner loop fills sites for each vector

		myPointer = &this->basis_vectors[i-simdata.M];	
		myPointer2 = &this->basis_vectors[i];

		cblas_dcopy(simdata.M, myPointer, 1, prev_basis, 1);
		cblas_dcopy(simdata.M, myPointer, 1, new_basis, 1);

		int z;
		int count = 0;
		nonzero = (double*)mkl_malloc((1 + count)*sizeof(double), 64);
		for (z = 0; z < simdata.M-1; ++z) {
			if (prev_basis[z] != 0) {
				mkl_realloc(nonzero, (1+count)*sizeof(double));
				nonzero[count] = z;
				count += 1;
			}
		}




		n_k = *std::max_element(nonzero, nonzero+count);

		if (n_k > 0) {
			cblas_dcopy(n_k-1, prev_basis, 1, new_basis, 1);
		}

		new_basis[n_k] = prev_basis[n_k] - 1;

		if (n_k + 1 <= simdata.M - 1) {
			new_basis[n_k + 1] = simdata.N - cblas_dasum(n_k+1, new_basis, 1);
		}

		if (n_k + 2 <= simdata.M - 1) {
			for (int z = n_k+2; z < simdata.M; ++z) {
				new_basis[z] = 0;
			}
		}

//		for (int j = 0; j < simdata.M; ++j) {
//			std::cout << new_basis[j];
//
//		}
//		std::cout << ";" << std::endl;


		cblas_dcopy(simdata.M, new_basis, 1, myPointer2, 1);

	}



	

	for (int i = 0; i < simdata.D; ++i) {
		for (int j = 0; j < simdata.M; ++j) {
			index = i*simdata.M + j;
			std::cout << basis_vectors[index];

		}
		std::cout << ";" << std::endl;
	}

};


Basis::~Basis(){
	mkl_free(basis_vectors);
	mkl_free(tags);
	mkl_free(tags_sorted);
};
