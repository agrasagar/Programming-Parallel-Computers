#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>

void mf(int ny, int nx, int hy, int hx, const float* in, float* out) {
    
	int xy;

   // medianFilter.reserve(hx*hy);
	#pragma omp parallel for schedule (static,1) shared(out) 
        for(xy = 0; xy < nx*ny; ++xy){   //for traversing through x axis
	
	int x = xy / ny;
        int y = xy % ny;
	std::vector<float> medianFilter;
	int size;
	float median;
	int window_y_min = std::max(y - hy, 0);
	int window_y_max = std::min(y + hy, (ny-1));
        int window_x_min = std::max(x - hx, 0);
	int window_x_max = std::min(x + hx, (nx-1));
        int window_y;
	int window_x;

	for(window_y = window_y_min; window_y<= window_y_max ; window_y++){ //traversing window y axis

            for(window_x = window_x_min ; window_x<= window_x_max ; window_x++){  //traversing window x axis

                  medianFilter.push_back(in[window_x + window_y*nx]); //filling up median filter vector

            }

        }

        size = medianFilter.size()/2;
        std::nth_element(medianFilter.begin(), medianFilter.begin() + (size), medianFilter.end());
        median = medianFilter.at(size);
        if ((medianFilter.size() % 2) == 0){
          std::nth_element(medianFilter.begin(), medianFilter.begin()+(size-1), medianFilter.begin()+size);
          median = (median + medianFilter[size-1])/2.0;

        }
        out[x + y*nx] = median;

      }

}

