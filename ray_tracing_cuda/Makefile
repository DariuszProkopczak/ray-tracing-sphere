all: ray_cuda

ray_cuda: ray_tracing.cu
	nvcc -o ray_cuda -arch=sm_37 -std=c++11 -lgomp -O3 ray_tracing.cu