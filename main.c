#include "mathlib.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define M 1.0        // black hole mass
#define Rs (2.0 * M) // Schwarzschild radius (rₛ) in geometrized units (G=c=1)

// define metric type
#define SCHWARZSCHILD_METRIC


/*
	this function doesn't take a Vector4 as argument as it won't be easy to calculate derivatives
	without modifying the vector itself. So, it directly uses Boyer-Lindquist coordinates, so
	adding a small value to a coordinate is direct.
*/
double calculate_metric_component(Vector4 x, int m, int n) {
	/*
		Schwarzschild line element (can directly derive the metric tensor):
		ds² = -(1 - rₛ/r)dt² + 1/(1 - rₛ/r)dr² + r²dθ² + r²sin²θdφ²
	*/
	double t, r, theta, phi;
	t = x.data[0];
	r = x.data[1];
	theta = x.data[2];
	phi = x.data[3];

	#ifdef SCHWARZSCHILD_METRIC
	if (m != n) return 0;
	if (n == 0) return -(1 - Rs / r);
	if (n == 1) return 1 / (1 - Rs / r);
	if (n == 2) return r * r;
	if (n == 3) return pow(r * sin(theta), 2);
	#elifdef MINKOWSKI_METRIC // for testing
	if (m != n) return 0;
	if (n == 0) return -1;
	else return 1;
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

#define STEP_SIZE 0.0001

State *get_next(State *s) {
	State *next = malloc(sizeof(State));
	next->k = malloc(sizeof(Vector4));
	next->x = malloc(sizeof(Vector4));

	Vector4 k, x;
	Vector4 a1, a2, a3, a4;
	Vector4 b1, b2, b3, b4;

	Metric *temp_g;
	Christoffel *temp_gamma;
	
	Vector4 *a = calculate_acceleration(s->gamma, s->k);
	
	// step 1
	for (int i = 0; i < 4; i++) {
		a1.data[i] = STEP_SIZE * a->data[i];
		b1.data[i] = STEP_SIZE * s->k->data[i];
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

	free(a);
	a = calculate_acceleration(temp_gamma, &k);

	for (int i = 0; i < 4; i++) {
		a2.data[i] = STEP_SIZE * a->data[i];
		b2.data[i] = STEP_SIZE * (s->k->data[i] + 0.5 * a1.data[i]);

		// prepare for step 3 by changing velocity vector based on a2
		k.data[i] = s->k->data[i] + a2.data[i] / 2.0;
	}

	// step 3
	/* 
		calculate new acceleration based on a2. position (therefore metric and christoffel) remain unchanged,
		since λ is not incremented
	
	*/
	free(a);
	a = calculate_acceleration(temp_gamma, &k);

	for (int i = 0; i < 4; i++) {
		a3.data[i] = STEP_SIZE * a->data[i];
		b3.data[i] = STEP_SIZE * (s->k->data[i] + 0.5 * a2.data[i]);

		// prepare for step 4 by incrementing position and velocity by a full step.
		k.data[i] = s->k->data[i] + a3.data[i];
		x.data[i] = s->x->data[i] + b3.data[i];
	}

	// step 4
	free(a);
	free(temp_g);
	free(temp_gamma);

	temp_g = calculate_metric(x);
	temp_gamma = calculate_christoffel(temp_g, &x);
	a = calculate_acceleration(temp_gamma, &k);

	for (int i = 0; i < 4; i++) {
		a4.data[i] = STEP_SIZE * a->data[i];
		b4.data[i] = STEP_SIZE + (s->k->data[i] * a3.data[i]);

		// calculate final velocity and position
		k.data[i] = s->k->data[i] + (a1.data[i] + 2 * a2.data[i] + 2 * a3.data[i] + a4.data[i]) / 6.0;
		x.data[i] = s->x->data[i] + (b1.data[i] + 2 * b2.data[i] + 2 * b3.data[i] + b4.data[i]) / 6.0;
	}

	// calculate final metric and christoffel symbols
	Metric *metric = calculate_metric(x);
	Christoffel *gamma = calculate_christoffel(metric, &x);

	// update next state
	next->g = metric;
	next->gamma = gamma;
	memcpy(next->k, &k, sizeof(Vector4));
	memcpy(next->x, &x, sizeof(Vector4));

	// free old state
	free(s->k);
	free(s->x);
	free(s->g);
	free(s->gamma);
	free(s);

	return next;
}

int main() {
	State *s = malloc(sizeof(State));
	s->x = malloc(sizeof(Vector4));
	s->k = malloc(sizeof(Vector4));

	// TODO: use photon null condition to calculate initial k^t
	s->x->data[0] = 0; // t
	s->x->data[1] = 10; // r
	s->x->data[2] = 1.57; // theta
	s->x->data[3] = 0; // phi
	
	s->k->data[0] = 1.25;
	s->k->data[1] = -1;
	s->k->data[2] = 0;
	s->k->data[3] = 0;
	
	s->g = calculate_metric(*(s->x));
	s->gamma = calculate_christoffel(s->g, s->x);

	State *next = get_next(s);

	printf(
		"STATE\nx = (%lf, %lf, %lf, %lf)\nk = (%lf, %lf, %lf, %lf)\ng = diag(%lf, %lf, %lf, %lf)\n",
		next->x->data[0], next->x->data[1], next->x->data[2], next->x->data[3],
		next->k->data[0], next->k->data[1], next->k->data[2], next->k->data[3],
		next->g->data[0][0], next->g->data[1][1], next->g->data[2][2], next->g->data[3][3]
	);
	return 0;
}
