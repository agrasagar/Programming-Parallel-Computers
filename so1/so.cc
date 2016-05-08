#include "so.h"
#include <math.h>
#include <algorithm>
#include <iostream>
#include <omp.h>
void psort(int n, data_t* data) {
    // FIXME: make this more efficient with parallelism
	int npt=omp_get_max_threads();
	int eachpt=n/npt;
	int tmp=npt;
	#pragma omp parallel for schedule(dynamic)
	for(int i=0;i<npt;i++){
		if(i==npt-1){
			std::sort(data+(npt-1)*eachpt,data+n);
		}
		else{
			std::sort(data+i*eachpt, data + (i+1)*eachpt);
		}
	}

	while(tmp!=1){
		#pragma omp parallel for schedule(dynamic)
		for(int i=0;i<tmp;i+=2){
			if((i+2)>=(tmp-1)){
				std::inplace_merge(data+i*eachpt,data+(i+1)*eachpt,data+n);
			}
			else{
				std::inplace_merge(data+i*eachpt,data+(i+1)*eachpt,data+(i+2)*eachpt);
			}
		}
		eachpt=eachpt*2;
		tmp=(tmp+1)*0.5;
	}

		
}
