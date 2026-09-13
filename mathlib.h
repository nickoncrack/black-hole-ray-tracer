#pragma once

typedef struct {
	double data[4][4];
} Metric;

typedef struct {
  	double data[4][4][4];
} Christoffel;

typedef struct {
  	double data[4];
} Vector4;

typedef struct {
  	Metric *g;
  	Christoffel *gamma;
  	Vector4 *x; // position
  	Vector4 *k; // velocity
} State;

/*
Metric *calculate_metric(State *s);
Metric *calculate_inverse_metric(Metric *g);

Christoffel *calculate_christoffel(Metric *g, Vector4 *pos);
State *get_next_state(State *s);
*/
