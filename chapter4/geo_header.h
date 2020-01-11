#pragma once

#include <iostream>
#include <vector>
#include <algorithm>

const double eps = 1e-8;

const double INF = 100;

int sign(const double &x); // Sign function. Avoid floating point error

class Point {
public:
	std::vector<double> val;

	// Constructors
	Point();
	Point(int d);
	Point(int d, std::vector<double> v);
	~Point();

	// Geometry operators
	Point operator + (const Point &p) const; // Addition
	Point operator - (const Point &p) const; // Sustraction
	double operator * (const Point &p) const; // inner product
	Point operator * (const double &k) const;
	double operator ^ (const Point &p) const; // inner product

	// Other functions
	int get_dim(); // Get the value of dimension
	void print_point(); // Print the point value

private:
	int dim; // Point dimension
};

class Polygon {
public:
	std::vector<Point> vertices;

	// Constructors
	Polygon();
	Polygon(int d);
	Polygon(int d, std::vector<Point> v);
	~Polygon();

	// Other functions
	int get_dim();

private:
	int dim; // Point dimension
};

class Line {
public:
	Point st, ed; // Two points on the line

	// Constructors
	Line();
	Line(int d, Point a, Point b);
	~Line();

	// Other functions
	int get_dim();

private:
	int dim; // Point dimension
};

Point line_intersection(Line &L1, Line &L2);

class HalfPlane {
public:
	// A HalfPlane is represented as: a[0] + ... + a[d-1] <= b
	Point a;
	double b;

	// Constructors
	HalfPlane();
	HalfPlane(int d, Point &p, double &b);
	~HalfPlane();

	// Other functions
	int get_dim();

private:
	int dim; // Point dimension
};

class LinearProgramming {
public:

	std::vector<HalfPlane> H; // a set of halfplanes
	Point O; // objective function
	
	// Constructors
	LinearProgramming();
	LinearProgramming(int d, std::vector<HalfPlane> &H, Point &O);
	LinearProgramming(int d);
	~LinearProgramming();

	std::vector<Point> solve();
	int get_dim();

private:
	int dim; // Point dimension
};

Line halfplane_boundary(HalfPlane &h);
