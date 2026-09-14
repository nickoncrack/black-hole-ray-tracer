#include "mathlib.h"

#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define a 1.0		 // black hole spin	
#define M 1.0        // black hole mass
#define Rs (2.0 * M) // Schwarzschild radius (rₛ) in geometrized units (G=c=1)

#define DIFF_STEP 1e-6

// define metric type
#define SCHWARZSCHILD_METRIC
// #define KERR_METRIC


double calculate_metric_component(Vector4 x, int m, int n) {
	/*
		Schwarzschild line element (can directly derive the metric tensor):
		ds² = -(1 - rₛ/r)dt² + 1/(1 - rₛ/r)dr² + r²dθ² + r²sin²θdφ²
	*/
	double t, r, θ, φ;
	t = x.data[0];
	r = x.data[1];
	θ = x.data[2];
	φ = x.data[3];

	#ifdef SCHWARZSCHILD_METRIC
	if (m != n) return 0;
	if (n == 0) return -(1 - Rs / r);
	if (n == 1) return 1 / (1 - Rs / r);
	if (n == 2) return r * r;
	if (n == 3) return pow(r * sin(θ), 2);
	#elifdef KERR_METRIC
	double ρ = pow(r, 2) + pow(a * cos(θ), 2);
	double ρ2 = pow(ρ, 2);

	if (m == n) {
		if (n == 0) { // dt²
			return -(1 - 2 * M * r / ρ2);
		} else if (n == 1) { // dr²
			return ρ2 / (pow(r, 2) - 2 * M * r + pow(a, 2));
		} else if (n == 2) { // dθ²
			return ρ2;
		} else if (n == 3) { // dφ²
			return (pow(r, 2) + pow(a, 2) + (2 * M * r * pow(a * sin(θ), 2) / ρ2)) * pow(sin(θ), 2);
		}
	} else {
		if ((m == 0 && n == 3) || (m == 3 && n == 0)) { // dtdφ
			return -1.0 * 4 * M * a * r * pow(sin(θ), 2) / ρ2;
		} else {
			return 0;
		}
	}
	#endif
}

Metric *calculate_metric(Vector4 x) {
	Metric *g = malloc(sizeof(Metric));
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			g->data[i][j] = calculate_metric_component(x, i, j);
		}
	}

  	return g;
}

/*
	the notation is horrible but i had to make this somehow readable

            gₘₙ(xᵞ + ε) - gₘₙ(xᵞ - ε)
	∂ᵧgₘₙ = -------------------------
	                   2ε
*/
double differentiate_metric_component(Metric *g, Vector4 pos, int m, int n, int gamma) {
	if (gamma == 0) return 0; // metric components don't change with respect to time on any metric
	
	#if defined(MINKOWSKI_METRIC) || defined(SCHWARZSCHILD_METRIC)
	if (m != n) return 0; // the schwarzschild metric doesn't have non-diagonal components
	#endif
	
	Vector4 pos_inc = pos;
	pos_inc.data[gamma] += DIFF_STEP;

	// use central difference to approximate the derivative which has O(ε²) truncation error
	double d1 = calculate_metric_component(pos_inc, m, n);
	
	pos_inc.data[gamma] -= 2 * DIFF_STEP;
	double d2 = calculate_metric_component(pos_inc, m, n);

	return (d1 - d2) / (2.0 * DIFF_STEP);
}


static double __det3(double m00, double m01, double m02, double m10, double m11,
                     double m12, double m20, double m21, double m22) {
  	return m00 * (m11 * m22 - m12 * m21) - m01 * (m10 * m22 - m12 * m20) +
           m02 * (m10 * m21 - m11 * m20);
}

Metric *calculate_inverse_metric(Metric *g) {
  	// any metric that will be used is a rank 2 tensor (matrix)
  	// Inverse matrix: A⁻¹ = 1/det(A) * Cᵀ, C: Cofactor matrix

	double A[4][4];
  	double C[4][4]; // cofactor matrix

  	memcpy(A, g->data, 16 * sizeof(double));

  	C[0][0] =  __det3(A[1][1], A[1][2], A[1][3], A[2][1], A[2][2], A[2][3], A[3][1], A[3][2], A[3][3]);
  	C[0][1] = -__det3(A[1][0], A[1][2], A[1][3], A[2][0], A[2][2], A[2][3], A[3][0], A[3][2], A[3][3]);
  	C[0][2] =  __det3(A[1][0], A[1][1], A[1][3], A[2][0], A[2][1], A[2][3], A[3][0], A[3][1], A[3][3]);
  	C[0][3] = -__det3(A[1][0], A[1][1], A[1][2], A[2][0], A[2][1], A[2][2], A[3][0], A[3][1], A[3][2]);

  	C[1][0] = -__det3(A[0][1], A[0][2], A[0][3], A[2][1], A[2][2], A[2][3], A[3][1], A[3][2], A[3][3]);
  	C[1][1] =  __det3(A[0][0], A[0][2], A[0][3], A[2][0], A[2][2], A[2][3], A[3][0], A[3][2], A[3][3]);
  	C[1][2] = -__det3(A[0][0], A[0][1], A[0][3], A[2][0], A[2][1], A[2][3], A[3][0], A[3][1], A[3][3]);
  	C[1][3] =  __det3(A[0][0], A[0][1], A[0][2], A[2][0], A[2][1], A[2][2], A[3][0], A[3][1], A[3][2]);

  	C[2][0] =  __det3(A[0][1], A[0][2], A[0][3], A[1][1], A[1][2], A[1][3], A[3][1], A[3][2], A[3][3]);
  	C[2][1] = -__det3(A[0][0], A[0][2], A[0][3], A[1][0], A[1][2], A[1][3], A[3][0], A[3][2], A[3][3]);
  	C[2][2] =  __det3(A[0][0], A[0][1], A[0][3], A[1][0], A[1][1], A[1][3], A[3][0], A[3][1], A[3][3]);
  	C[2][3] = -__det3(A[0][0], A[0][1], A[0][2], A[1][0], A[1][1], A[1][2], A[3][0], A[3][1], A[3][2]);
  	
	C[3][0] = -__det3(A[0][1], A[0][2], A[0][3], A[1][1], A[1][2], A[1][3], A[2][1], A[2][2], A[2][3]);
  	C[3][1] =  __det3(A[0][0], A[0][2], A[0][3], A[1][0], A[1][2], A[1][3], A[2][0], A[2][2], A[2][3]);
  	C[3][2] = -__det3(A[0][0], A[0][1], A[0][3], A[1][0], A[1][1], A[1][3], A[2][0], A[2][1], A[2][3]);
  	C[3][3] =  __det3(A[0][0], A[0][1], A[0][2], A[1][0], A[1][1], A[1][2], A[2][0], A[2][1], A[2][2]);

  	/*
		Calculate determinant using Laplace expansion on the first row
		det A = sum from j=0 to 3 of A₀ⱼC₀ⱼ
  	*/
  	double det = A[0][0] * C[0][0] + A[0][1] * C[0][1] + A[0][2] * C[0][2] + A[0][3] * C[0][3];
  	
  	if (fabs(det) < 1e-12) {
  		// metric is not invertible (unlikely scenario since all metric tensors are invertible)
  		printf("Metric is not invertible\n");
  		return NULL;
  	}

  	double det_inv = 1.0 / det;
	Metric *g_inv = malloc(sizeof(Metric));

	// scale the transpose of the cofactor matrix by 1/det(A)
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			g_inv->data[i][j] = C[i][j] * det_inv;
		}
	}
	
	return g_inv;
}


double __calculate_christoffel_component(Metric *g, Metric *g_inv, Vector4 *pos, int i, int j, int k) {
	double ret = 0;

	for (int m = 0; m < 4; m++) {
		// if (fabs(g_inv->data[k][m]) < 1e-12) continue;

		double d1 = differentiate_metric_component(g, *pos, m, j, i); // ∂ᵢgₘⱼ 
		double d2 = differentiate_metric_component(g, *pos, m, i, j); // ∂ⱼgₘᵢ
		double d3 = differentiate_metric_component(g, *pos, i, j, m); // ∂ₘgᵢⱼ

		ret += g_inv->data[k][m] * (d1 + d2 - d3);
	}

	return 0.5 * ret;
}

/*
	Γᵏᵢⱼ = 1/2 * gᵐᵏ(∂ᵢgₘⱼ + ∂ⱼgₘᵢ - ∂ₘgᵢⱼ)
	
	metric tensors are symmetric, so christoffel symbols are also symmetric:
	Γᵏᵢⱼ = Γᵏⱼᵢ
*/
Christoffel *calculate_christoffel(Metric *g, Vector4 *pos) {
	Christoffel *gamma = malloc(sizeof(Christoffel));
	Metric *g_inv = calculate_inverse_metric(g);

	for (int k = 0; k < 4; k++) {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				gamma->data[k][i][j] = __calculate_christoffel_component(
					g, g_inv, pos, i, j, k
				);
			}
		}
	}

	free(g_inv);
	return gamma;
}

// a = -Γᵏᵢⱼ kⁱkʲ
double calculate_acceleration_component(Christoffel *gamma, Vector4 *v, int k) {
	double ret = 0.0;

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			ret += gamma->data[k][i][j] * v->data[i] * v->data[j];
		}
	}

	return ret;
}

Vector4 *calculate_acceleration(Christoffel *gamma, Vector4 *k) {
	Vector4 *ret = malloc(sizeof(Vector4));

	for (int i = 0; i < 4; i++) ret->data[i] = -1.0 * calculate_acceleration_component(gamma, k, i);

	return ret;
}


#define BASE_STEP_SIZE 0.01

// using a step size less than approximately 0.0001 causes the output numbers to have a huge error
double get_step_size(double r) {
	double factor = (r - 2.0 * M) / (8.0 * M);

	if (factor < 0.1) factor = 0.1;
	if (factor > 1.0) factor = 1;

	return BASE_STEP_SIZE * factor;
}

State *get_next(State *s, double Δλ) {
	State *next = malloc(sizeof(State));
	next->k = malloc(sizeof(Vector4));
	next->x = malloc(sizeof(Vector4));

	Vector4 k, x;
	Vector4 a1, a2, a3, a4;
	Vector4 b1, b2, b3, b4;

	Metric *temp_g;
	Christoffel *temp_gamma;
	
	Vector4 *ac = calculate_acceleration(s->gamma, s->k);
	
	// step 1
	for (int i = 0; i < 4; i++) {
		a1.data[i] = Δλ * ac->data[i];
		b1.data[i] = Δλ * s->k->data[i];
	}

	// step 2
	// increment velocity and position vectors by a half time step
	for (int i = 0; i < 4; i++) {
		k.data[i] = s->k->data[i] + a1.data[i] / 2.0;
		x.data[i] = s->x->data[i] + b1.data[i] / 2.0;
	}
 
	// calculate new metric and christoffel symbols
	temp_g = calculate_metric(x);
	temp_gamma = calculate_christoffel(temp_g, &x);

	free(ac);
	ac = calculate_acceleration(temp_gamma, &k);

	for (int i = 0; i < 4; i++) {
		a2.data[i] = Δλ * ac->data[i];
		b2.data[i] = Δλ * (s->k->data[i] + 0.5 * a1.data[i]);

		// prepare for step 3 by changing velocity vector based on a2
		k.data[i] = s->k->data[i] + a2.data[i] / 2.0;
		x.data[i] = s->x->data[i] + b2.data[i] / 2.0;
	}

	free(ac);
	free(temp_g);
	free(temp_gamma);
	temp_g = calculate_metric(x);
	temp_gamma = calculate_christoffel(temp_g, &x);
	ac = calculate_acceleration(temp_gamma, &k);

	// step 3
	for (int i = 0; i < 4; i++) {
		a3.data[i] = Δλ * ac->data[i];
		b3.data[i] = Δλ * (s->k->data[i] + 0.5 * a2.data[i]);

		// prepare for step 4 by incrementing position and velocity by a full step.
		k.data[i] = s->k->data[i] + a3.data[i];
		x.data[i] = s->x->data[i] + b3.data[i];
	}

	// step 4
	free(ac);
	free(temp_g);
	free(temp_gamma);

	temp_g = calculate_metric(x);
	temp_gamma = calculate_christoffel(temp_g, &x);
	ac = calculate_acceleration(temp_gamma, &k);

	for (int i = 0; i < 4; i++) {
		a4.data[i] = Δλ * ac->data[i];
		b4.data[i] = Δλ * (s->k->data[i] + a3.data[i]);

		// calculate final velocity and position
		k.data[i] = s->k->data[i] + (a1.data[i] + 2 * a2.data[i] + 2 * a3.data[i] + a4.data[i]) / 6.0;
		x.data[i] = s->x->data[i] + (b1.data[i] + 2 * b2.data[i] + 2 * b3.data[i] + b4.data[i]) / 6.0;
	}

	free(ac);
	free(temp_g);
	free(temp_gamma);

	// calculate final metric and christoffel symbols
	Metric *metric = calculate_metric(x);
	Christoffel *gamma = calculate_christoffel(metric, &x);

	// update next state
	next->g = metric;
	next->gamma = gamma;
	memcpy(next->k, &k, sizeof(Vector4));
	memcpy(next->x, &x, sizeof(Vector4));

	return next;
}

void print_state(State *s) {
	printf(
		"STATE\nx = (%lf, %lf, %lf, %lf)\nk = (%lf, %lf, %lf, %lf)\ng = diag(%lf, %lf, %lf, %lf)\n",
		s->x->data[0], s->x->data[1], s->x->data[2], s->x->data[3],
		s->k->data[0], s->k->data[1], s->k->data[2], s->k->data[3],
		s->g->data[0][0], s->g->data[1][1], s->g->data[2][2], s->g->data[3][3]
	);

	return;
}

void free_state(State *s) {
	free(s->x);
	free(s->k);
	free(s->g);
	free(s->gamma);
	free(s);

	return;
}

enum photon_state {
	P_TERMINATE,
	P_ESCAPE,
	P_HIT_DISK
};

unsigned long long steps = 0;
unsigned long long terminated_photons = 0;

#define R_IN (3.0 * Rs) // ISCO
#define R_OUT (15.0 * M)

enum photon_state get_final_state(Vector4 *x, Vector4 *k, Vector4 **final_x) {
	State *s = malloc(sizeof(State));
	State *new_state;
	s->x = x;
	s->k = k;
	s->g = calculate_metric(*x);
	s->gamma = calculate_christoffel(s->g, x);

	int step_count = 0;
	while (step_count < 10000) {
		double r = s->x->data[1];
		new_state = get_next(s, get_step_size(r));
		r = new_state->x->data[1];

		// check if photon crossed the accretion disk
		double theta_old = s->x->data[2] - M_PI / 2.0;
		double theta_new = new_state->x->data[2] - M_PI / 2.0;

		free_state(s); // free old state
		s = new_state;
		if (theta_old * theta_new < 0) {
			if (r >= R_IN && r <= R_OUT) {
				free_state(s);
				return P_HIT_DISK;
			}
		}

		// slightly increase termination radius to eliminate noise
		if (s->x->data[1] <= 1.025 * Rs || isnan(r)) {
			free_state(s);

			#pragma omp atomic
			terminated_photons++;

			return P_TERMINATE;
		} else if (s->x->data[1] > 50 * M) {
			goto escape;
		}

		#pragma omp atomic
		steps++;

		step_count++;
	}

escape:
	memcpy(*final_x, s->x, sizeof(Vector4)); // for color data
	free_state(s);

	return P_ESCAPE;
}

#define SCREEN_WIDTH	400
#define SCREEN_HEIGHT   400
#define CAMERA_FOV	(4.0 * M_PI / 9.0) // 60deg
#define K2(n) pow(k->data[n], 2)

Vector4 *calculate_velocity_from_pixel_coords(int i, int j, Vector4 *camera_pos) {
	// normalize coordinates
	double u = (i + 0.5) / SCREEN_WIDTH - 0.5;
	double v = (j + 0.5) / SCREEN_HEIGHT - 0.5;

	double fov_tan = tan(CAMERA_FOV / 2.0);
	double x_screen = u * fov_tan;
	double y_screen = -v * fov_tan * SCREEN_HEIGHT / SCREEN_WIDTH;

	double g_tt = calculate_metric_component(*camera_pos, 0, 0);
	double g_rr = calculate_metric_component(*camera_pos, 1, 1);
	double g_θθ = calculate_metric_component(*camera_pos, 2, 2);
	double g_φφ = calculate_metric_component(*camera_pos, 3, 3);

	Vector4 *k = malloc(sizeof(Vector4));
	k->data[1] = - 1.0 / sqrt(g_rr);
	k->data[2] = y_screen / sqrt(g_θθ);
	k->data[3] = x_screen / sqrt(g_φφ);

	// use null condition to calculate k^t
	k->data[0] = sqrt(-1/g_tt * (g_rr * K2(1) + g_θθ * K2(2) + g_φφ * K2(3)));

	return k;
}

// generate a checkerboard pattern as background
void get_bg_color(double theta, double phi, Color *color) {
	phi = fmod(phi, 2.0 * M_PI);
	if (phi < 0) phi += 2.0 * M_PI;

	double tile_size = M_PI / 12.0;
	int theta_idx = (int) (theta / tile_size);
	int phi_idx = (int) (phi / tile_size);

	if ((theta_idx + phi_idx) % 2 == 0) {
		// light gray
		color->r = 15;
		color->g = 20;
		color->b = 30;
	} else {
		// light blue
		color->r = 30;
		color->g = 42;
		color->b = 60;
	}

	return;
}

void save_ppm(const char *filename, unsigned char *buff) {
	FILE *fp = fopen(filename, "wb");
	if (!fp) {
		perror("Failed to open file\n");
		return;
	}

	fprintf(fp, "P6\n%d %d\n255\n", SCREEN_WIDTH, SCREEN_HEIGHT);
	fwrite(buff, 3, SCREEN_WIDTH * SCREEN_HEIGHT, fp);
	fclose(fp);

	printf("Saved image to %s\n", filename);
	return;
}

double get_time() {
	struct timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	return (double) ts.tv_sec + (double) ts.tv_nsec * 1e-9;
}

int main() {
	Vector4 camera_pos;
	camera_pos.data[0] = 0;			// t
	camera_pos.data[1] = 50 * M;	// r
	camera_pos.data[2] = (5.0 * M_PI / 12.0);  // θ =75deg, slightly looking down at the accretion disk
	camera_pos.data[3] = 0;			// φ

	unsigned char *buff = malloc(SCREEN_WIDTH * SCREEN_HEIGHT * 3);

	double start, stop;
	start = get_time();

	#pragma omp parallel for schedule(dynamic)
	for (int y = 0; y < SCREEN_HEIGHT; y++) {
		for (int x = 0; x < SCREEN_WIDTH; x++) {
			Vector4 *photon_x = malloc(sizeof(Vector4));
			Vector4 *final_x = malloc(sizeof(Vector4));
			memcpy(photon_x, &camera_pos, sizeof(Vector4));

			Vector4 *k = calculate_velocity_from_pixel_coords(x, y, photon_x);
			enum photon_state final = get_final_state(photon_x, k, &final_x);

			int idx = (y * SCREEN_WIDTH + x) * 3;
			if (final == P_TERMINATE) {
				buff[idx + 0] = 0;
				buff[idx + 1] = 0;
				buff[idx + 2] = 0;
			} else if (final == P_ESCAPE) {
				Color color;
				get_bg_color(final_x->data[2], final_x->data[3], &color);

				buff[idx + 0] = color.r;
				buff[idx + 1] = color.g;
				buff[idx + 2] = color.b;
			} else { // accretion disk
				buff[idx + 0] = 255;
				buff[idx + 1] = 138;
				buff[idx + 2] = 61;
			}

			free(final_x);
		}
	}

	stop = get_time();

	printf("steps: %llu\nterminated photons: %llu\ntime elapsed: %lf s\n", steps, terminated_photons, stop - start);

	save_ppm("black_hole.ppm", buff);
	free(buff);

	return 0;
}
