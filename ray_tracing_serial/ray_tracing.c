#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#define M_PI 3.14159265358979323846264338327

typedef struct vector {
	double x;
	double y;
	double z;
} Vector;

Vector create_vector(double x, double y, double z)
{
	Vector vect;
	vect.x = x;
	vect.y = y;
	vect.z = z;
	return vect;
}

double dot_product(Vector * in_vect1, Vector * in_vect2)
{
	double out = (in_vect1->x * in_vect2->x);
	out = out + (in_vect1->y * in_vect2->y);
	out = out + (in_vect1->z * in_vect2->z);

	return out;
}

double norm(Vector * in_vect)
{
	double out = (in_vect->x * in_vect->x) + (in_vect->y * in_vect->y) + (in_vect->z * in_vect->z);
	out = sqrt(out);

	return out;
}

void scalar_multiply(Vector * in_vect, Vector * out_vect, double scalar)
{
	out_vect->x = in_vect->x * scalar;
	out_vect->y = in_vect->y * scalar;
	out_vect->z = in_vect->z * scalar;
}

void scalar_divide(Vector * in_vect, Vector * out_vect, double scalar)
{
	out_vect->x = in_vect->x / scalar;
	out_vect->y = in_vect->y / scalar;
	out_vect->z = in_vect->z / scalar;
}

void vector_add(Vector * in_vect1, Vector * in_vect2, Vector * out_vect)
{
	out_vect->x = in_vect1->x + in_vect2->x;
	out_vect->y = in_vect1->y + in_vect2->y;
	out_vect->z = in_vect1->z + in_vect2->z;
}

void vector_subtract(Vector * in_vect1, Vector * in_vect2, Vector * out_vect)
{
	out_vect->x = in_vect1->x - in_vect2->x;
	out_vect->y = in_vect1->y - in_vect2->y;
	out_vect->z = in_vect1->z - in_vect2->z;
}

void direction_sample(Vector * vect)
{

	double phi = ((double) rand() / RAND_MAX) * 2 * M_PI;
	double cos_theta = (((double) rand() / RAND_MAX) * 2) - 1;
	double sin_theta = sqrt(1 - (cos_theta * cos_theta));
	vect->x = sin_theta * cos(phi);
	vect->y = sin_theta * sin(phi);
	vect->z = cos_theta;
}

double max(double a, double b)
{
	if (a >= b)
		return a;
	else
		return b;
}

void write_animation(double ** animation_data, int number_of_frames, int num_grid_points, char * file_name)
{

	FILE *fp = fopen(file_name, "w");

	for (int frame_number=0; frame_number<number_of_frames; frame_number++)
	{
		for (int row_idx=(num_grid_points-1); row_idx>=0; row_idx--)
		{
			for(int col_idx=0; col_idx<num_grid_points; col_idx++)
			{
			 	fprintf(fp, "%f ", animation_data[frame_number][col_idx + (num_grid_points * row_idx)]);
			}

			fprintf(fp, "\n");
		}
	}

	fclose(fp);

}

void write_G(double * G, int num_grid_points, char * file_name)
{
	// writes data to file
	FILE *fp = fopen(file_name, "w");

	for (int row_idx=(num_grid_points-1); row_idx>=0; row_idx--)
	{
		for(int col_idx=0; col_idx<num_grid_points; col_idx++)
		 	fprintf(fp, "%f ", G[col_idx + (num_grid_points * row_idx)]);
		fprintf(fp, "\n");
	}

	fclose(fp);

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

	// animation parameters
	int create_animation = 0;											// if false, only one frame is saved
	int number_of_frames = 40;											// number of animation frames
	Vector L_v = create_vector(0, 40, 5);								// velocity of light source

	// create animation data structure
	double ** G_animation;
	if (create_animation)
	{
		G_animation = malloc(sizeof(double *) * number_of_frames);
		for (int i=0; i<number_of_frames; i++)
		{
			G_animation[i] = malloc(sizeof(double) * num_grid_points * num_grid_points);
			for (int j=0; j<(num_grid_points * num_grid_points); j++)
				G_animation[i][j] = 0;
		}
	}
	else
		number_of_frames = 1;

	// random seed for generating light ray directions
	srand(3);

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
	int row, col;										// grid points of G representing W

	// create grid and set all values to 0
	double * G = malloc(sizeof(double) * num_grid_points * num_grid_points);
	for (int i=0; i<(num_grid_points*num_grid_points); i++)
		G[i] = 0;

	// loop through frames if performing animation, otherwise just run one frame
	for (int frame_idx = 0; frame_idx < number_of_frames; frame_idx++)
	{
		//main algorithm
		for (int i=0; i<number_of_rays; i++)
		{

			temp = -1;	// initialize to value less than 0 so loop is entered

			while ((W.x > W_max) || (W.x < (-1*W_max)) || (W.z > W_max)  || (W.z < (-1*W_max)) || (temp <= 0))
			{
				direction_sample(&V);
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

			b = max(0, dot_product(&S, &N));

			row = round(num_grid_points * (W.x + W_max) / (2*W_max));
			col = round(num_grid_points * (W.z + W_max) / (2*W_max));

			if (create_animation)
				G_animation[frame_idx][col + (num_grid_points * row)] = G_animation[frame_idx][col + (num_grid_points * row)] + b;
			else
				G[col + (num_grid_points * row)] = G[col + (num_grid_points * row)] + b;

		}

		if (create_animation)
			// move light source
			vector_add(&L, &L_v, &L);
	}

	// print time
	end_time = omp_get_wtime();
	printf("Total Run Time: %f secs\n", end_time - start_time);

	// write to file
	if (create_animation)
	{
		// write animation
		printf("Writing animation\n");
		write_animation(G_animation, number_of_frames, num_grid_points, "animation.ascii");
		for (int i=0; i<number_of_frames; i++)
			free(G_animation[i]);
		free(G_animation);
	}
	else
	{
		// write G
		printf("Writing G\n");
		write_G(G, num_grid_points, "plot.bin");
		free(G);
	}

	return 0;
	
}
