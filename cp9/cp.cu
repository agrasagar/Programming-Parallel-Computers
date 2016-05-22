#include "cp.h"
#include <iostream>     // std::cout
#include <algorithm>    // std::for_each
#include <vector>       // std::vector
#include <cmath>
#define TILE 32

__global__ void matrixMul(double *Md,double *Mtd, float *Rd, int ny, int nx){
	
	int x = blockIdx.x*TILE + threadIdx.x;
	int y = blockIdx.y*TILE + threadIdx.y;
	__shared__ float tile_A[TILE][TILE];
	__shared__ float tile_B[TILE][TILE];
	int i;
	float sum = 0.0;
	for(i = 0; i< ((nx-1)/TILE)+1; ++i ){
		
		//Put tile_A on shared memory
		if((i*TILE + threadIdx.x )< nx && y < ny ){
			tile_A[threadIdx.y][threadIdx.x] = Md[y*nx + i*TILE + threadIdx.x];
		}
		else{
			tile_A[threadIdx.y][threadIdx.x] = 0;
		}
		//Put tile_B on shared memory
		if((i*TILE + threadIdx.y)< nx && x < ny ){
			tile_B[threadIdx.y][threadIdx.x] = Mtd[(i*TILE + threadIdx.y)*ny + x];
		}
		else{
			tile_B[threadIdx.y][threadIdx.x] = 0;
		}
		__syncthreads();
		
		for(int j = 0; j < TILE; ++j){
			sum += tile_A[threadIdx.y][j] * tile_B[j][threadIdx.x];
		}
		__syncthreads();
	}
	
	if(x >= ny || y >= ny || x<y){
		return;
	}
	
	Rd[((blockIdx.y*blockDim.y + threadIdx.y)*ny) + (blockIdx.x*blockDim.x) + threadIdx.x] = sum;
}

void correlate(int ny, int nx, const float* data, float* result) {
		
	double* dataVec;
	double* dataVecT;
	cudaMallocHost( (void **) &dataVec, sizeof(double)*nx*ny);
	cudaMallocHost( (void **) &dataVecT, sizeof(double)*ny*nx);
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
	//Matrix transposel
	for( y = 0 ; y <ny ; ++y){
		for(i = 0 ; i < nx; ++i){
			dataVecT[i*ny + y] = dataVec[y*nx + i];
		}
	} 
	//Matrix multiplication
	//matrix nx*ny
	//matrixT ny*nx
	//after multiplication product matrix would be of size ny*ny
	dim3 dimBlock(TILE,TILE);
	dim3 dimGrid(std::ceil((double)ny/dimBlock.x), std::ceil((double)ny/(double)dimBlock.y));
	int sizeM = nx*ny*sizeof(double);
	int sizeR = ny*ny*sizeof(float);
	double* Md = 0;
	double* Mtd = 0; 
	float* Rd = 0;
	//Allocate Md and Pd on device
	cudaMalloc((void **)&Md, sizeM);
	cudaMemcpy(Md, dataVec, sizeM, cudaMemcpyHostToDevice);
	cudaMalloc((void **)&Mtd, sizeM);
        cudaMemcpy(Mtd, dataVecT, sizeM, cudaMemcpyHostToDevice);

	cudaMalloc((void **)&Rd, sizeR);
	matrixMul<<<dimGrid, dimBlock>>>(Md,Mtd,Rd,ny, nx);
	cudaMemcpy(result, Rd, sizeR, cudaMemcpyDeviceToHost);

	cudaFree(Md);
	cudaFree(Mtd);
	cudaFree(Rd);
	cudaFree(dataVec);
	cudaFree(dataVecT);	
}
