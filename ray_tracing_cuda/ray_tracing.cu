#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <float.h>
#include <math.h>
#include <omp.h>
#define THREADS_PER_BLOCK 512

typedef struct vector {
	double x;
	double y;
	double z;
} Vector;

// A 63−bit LCG
// Returns a double precision value from a uniform distribution
// between 0.0 and 1.0 using a caller −owned state variable . 
__device__ double LCG_random_double(uint64_t * seed)
{
	const uint64_t m = 9223372036854775808ULL; // 2ˆ63
	const uint64_t a = 2806196910506780709ULL; 
	const uint64_t c = 1ULL;
	 
	*seed = (a * (*seed) + c) %m;
	return (double) (*seed) / (double) m;
}

__device__ Vector create_vector(double x, double y, double z)
{
	Vector vect;
	vect.x = x;
	vect.y = y;
	vect.z = z;
	return vect;
}

__device__ double dot_product(Vector * in_vect1, Vector * in_vect2)
{
	double out = (in_vect1->x * in_vect2->x);
	out = out + (in_vect1->y * in_vect2->y);
	out = out + (in_vect1->z * in_vect2->z);

	return out;
}

__device__ double norm(Vector * in_vect)
{
	double out = (in_vect->x * in_vect->x) + (in_vect->y * in_vect->y) + (in_vect->z * in_vect->z);
	out = sqrt(out);

	return out;
}

__device__ void scalar_multiply(Vector * in_vect, Vector * out_vect, double scalar)
{
	out_vect->x = in_vect->x * scalar;
	out_vect->y = in_vect->y * scalar;
	out_vect->z = in_vect->z * scalar;
}

__device__ void scalar_divide(Vector * in_vect, Vector * out_vect, double scalar)
{
	out_vect->x = in_vect->x / scalar;
	out_vect->y = in_vect->y / scalar;
	out_vect->z = in_vect->z / scalar;
}

__device__ void vector_add(Vector * in_vect1, Vector * in_vect2, Vector * out_vect)
{
	out_vect->x = in_vect1->x + in_vect2->x;
	out_vect->y = in_vect1->y + in_vect2->y;
	out_vect->z = in_vect1->z + in_vect2->z;
}

__device__ void vector_subtract(Vector * in_vect1, Vector * in_vect2, Vector * out_vect)
{
	out_vect->x = in_vect1->x - in_vect2->x;
	out_vect->y = in_vect1->y - in_vect2->y;
	out_vect->z = in_vect1->z - in_vect2->z;
}

__device__ void direction_sample(Vector * vect, uint64_t * seed)
{

	double rand1 = LCG_random_double(seed);
	double rand2 = LCG_random_double(seed);

	double phi = rand1 * 2 * M_PI;
	double cos_theta = (rand2 * 2) - 1;
	double sin_theta = sqrt(1 - (cos_theta * cos_theta));
	vect->x = sin_theta * cos(phi);
	vect->y = sin_theta * sin(phi);
	vect->z = cos_theta;
}

void write_G(float * G, int num_grid_points)
{
	// writes data to file
	FILE *fp = fopen("plot.bin", "w");

	for (int row_idx=(num_grid_points-1); row_idx>=0; row_idx--)
	{
		for(int col_idx=0; col_idx<num_grid_points; col_idx++)
		 	fprintf(fp, "%f ", G[col_idx + (num_grid_points * row_idx)]);
		fprintf(fp, "\n");
	}

	fclose(fp);

}

__global__ void ray_tracing(float * G, int num_grid_points, int W_max, int W_y, int R, int C_y, int L_x, int L_y, int L_z)
{ 

	// initialize vectors and other variables
	Vector V = create_vector(0, 0, 0);					// light ray direction
	Vector C = create_vector(0, C_y, 0);				// center of sphere
	Vector W = create_vector(0, W_y, 0);				// point where light ray intersects window
	Vector I = create_vector(0, 0, 0);					// point where light ray intersects sphere
	Vector N = create_vector(0, 0, 0);					// unit normal vector at I
	Vector L = create_vector(L_x, L_y, L_z);			// location of light source
	Vector S = create_vector(0, 0, 0);					// direction of light source at I
	Vector temp_vector;									
	double temp;										// represents (V*C)^2 + R^2 - C*C		
	double t;											// t
	double b;											// brightness at I
	int row, col;			

	uint64_t seed = (uint64_t) ((threadIdx.x + blockIdx.x*blockDim.x) * 5356812);

	temp = -1;	// initialize to value less than 0 so loop is entered

	while ((W.x > W_max) || (W.x < (-1*W_max)) || (W.z > W_max) || (W.z < (-1*W_max)) || (temp <= 0))
	{
		direction_sample(&V, &seed);
		scalar_multiply(&V, &W, W_y / V.y);
		temp = dot_product(&V, &C);
		temp = temp * temp;
		temp = temp + (R*R);
		temp = temp - dot_product(&C, &C);
	}

	t = dot_product(&V, &C) - sqrt(temp);
	scalar_multiply(&V, &I, t);

	vector_subtract(&I, &C, &temp_vector);
	scalar_divide(&temp_vector, &N, norm(&temp_vector));

	vector_subtract(&L, &I, &temp_vector);
	scalar_divide(&temp_vector, &S, norm(&temp_vector));

	b = dot_product(&S, &N);
	if (b<0)
		b = 0;

	row = round(num_grid_points * (W.x + W_max) / (2*W_max));
	col = round(num_grid_points * (W.z + W_max) / (2*W_max));

	atomicAdd(&G[col + (num_grid_points * row)], (float) b);

}

int main(int argc, char * argv[])
{

	if (argc != 3)
	{
		printf("2 arguments needed: number of rays, and number of grid points. Exiting.\n");
		return 0;
	}
	int number_of_rays = atoi(argv[1]);
	int num_grid_points = atoi(argv[2]);

	// set timing variables
	double start_time = omp_get_wtime();
	double end_time;

	// set hard-coded parameters
	double W_max = 100;					// width of window
	double W_y = 100;					// window y position
	double R = 60;						// radius of sphere
	double C_y = 200;					// sphere y position 
	double L_x = 500;					// light source x position
	double L_y = -100;					// light source y position
	double L_z = -400;					// light source z position

	// create grid and set all values to 0
	float * G = (float *) malloc(sizeof(float) * num_grid_points * num_grid_points);
	for (int i=0; i<(num_grid_points*num_grid_points); i++)
		G[i] = 0;
	float * dev_G;

	// set device memory
	cudaMalloc((void **) &dev_G, sizeof(float) * num_grid_points * num_grid_points);
	printf("cudaMalloc status: %s\n", cudaGetErrorString(cudaGetLastError()));
	cudaMemcpy(dev_G, G, sizeof(float) * num_grid_points * num_grid_points, cudaMemcpyHostToDevice);
	printf("cudaMemcpy status: %s\n", cudaGetErrorString(cudaGetLastError()));

  	/* launch number_of_rays threads organized in blocks of size THREADS_PER_BLOCK */
  	ray_tracing<<<number_of_rays/THREADS_PER_BLOCK,THREADS_PER_BLOCK>>>(dev_G, num_grid_points, W_max, W_y, R, C_y, L_x, L_y, L_z);
 	
	printf("kernel: %s\n", cudaGetErrorString(cudaGetLastError())); 

	cudaDeviceSynchronize();

	printf("synchronize status: %s\n", cudaGetErrorString(cudaGetLastError())); 

  	// copy back memory
  	cudaMemcpy(G, dev_G, sizeof(float) * num_grid_points * num_grid_points, cudaMemcpyDeviceToHost);
  	printf("cudaMemcpy status: %s\n", cudaGetErrorString(cudaGetLastError())); 

	// print time
	end_time = omp_get_wtime();
	printf("Total Run Time: %f secs\n", end_time - start_time);

	// write G
	printf("Writing G\n");
	write_G(G, num_grid_points);
	free(G);

	return 0;

}

