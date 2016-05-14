#include "cp.h"
#include <iostream>     // std::cout
#include <algorithm>    // std::for_each
#include <vector>       // std::vector
#include <cmath>
__global__ void matrixMul(double *Md, float *Rd, int ny, int nx){
	
	int x = blockIdx.x*blockDim.x + threadIdx.x;
	int y = blockIdx.y*blockDim.y + threadIdx.y;
	double sum = 0;	
	if(x >= ny || y >= ny || x<y){
		return;
	}
	for(int k= 0; k<nx ; k++){
		sum += Md[k + y*nx] * Md[k + x*nx];
	}
	Rd[x + y*ny] = sum;
}

void correlate(int ny, int nx, const float* data, float* result) {
		
	double* dataVec;
	cudaMallocHost( (void **) &dataVec, sizeof(double)*nx*ny);
	int y, rowStart, rowEnd;
	int i;
	double meanOfEachRow, meanOfSquareRoot,num;
	std::vector<double> v(nx), vzeroSquaremean(nx);
	for(y = 0; y< ny ; ++y){        //traversing through y axis
		rowStart = y*nx;
		rowEnd = rowStart + nx;
		meanOfEachRow = std::accumulate(data+rowStart, data+rowEnd, 0.0)/nx;
		std::transform(data+rowStart, data+rowEnd, v.begin(), [&meanOfEachRow](double val){ return (val - meanOfEachRow); });
		std::transform(v.begin(), v.end(), vzeroSquaremean.begin(), [](double val){ return std::pow(val, 2); });
		meanOfSquareRoot = std::sqrt(std::accumulate(vzeroSquaremean.begin(), vzeroSquaremean.end(), 0.0));
		for(i = 0; i < nx; ++i) {
			num = v[i]/meanOfSquareRoot;
			dataVec[i+y*nx] = num;
		}
	}
	//Matrix multiplication
	//matrix nx*ny
	//matrixT ny*nx
	//after multiplication product matrix would be of size ny*ny
	dim3 dimBlock(8,8);
	dim3 dimGrid(std::ceil((double)ny/dimBlock.x), std::ceil((double)ny/(double)dimBlock.y));
	int sizeM = nx*ny*sizeof(double);
	int sizeR = ny*ny*sizeof(float);
	double* Md = 0; 
	float* Rd = 0;
	//Allocate Md and Pd on device
	cudaMalloc((void **)&Md, sizeM);
	cudaMemcpy(Md, dataVec, sizeM, cudaMemcpyHostToDevice);
	cudaMalloc((void **)&Rd, sizeR);
	matrixMul<<<dimGrid, dimBlock>>>(Md, Rd,ny, nx);
	cudaMemcpy(result, Rd, sizeR, cudaMemcpyDeviceToHost);

	cudaFree(Md);
	cudaFree(Rd);
	cudaFree(dataVec);
	
}
