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

typedef struct {
	unsigned char r, g, b;
} Color;
