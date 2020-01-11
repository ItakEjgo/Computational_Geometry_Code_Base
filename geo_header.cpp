#include "geo_header.h"
#include <assert.h>


#pragma region Point

Point::Point() {}

Point::Point(int d) {
	dim = d;
	val.resize(d, 0);
}

Point::Point(int d, std::vector<double> v) {
	dim = d;
	val.resize(0);
	for (double Xi : v) 
		val.push_back(Xi);
}

Point::~Point() {}

// Get point dimension, dim is private
int Point::get_dim() {
	return dim;
}

// Point print, show dim and val
void Point::print_point() {
	std::cout << "dim = " << dim << std::endl;
	std::cout << "val =";
	for (double Xi : val) 
		std::cout << " " << Xi;
	std::cout << std::endl;
}

// Point add operator, left + right
Point Point::operator + (const Point &p) const {
	assert(dim == p.dim);
	std::vector<double> v = {};
	for (int i = 0; i < dim; i++) 
		v.push_back(val[i] + p.val[i]);
	return Point(dim, v);
}

// Point substract operator, left - right
Point Point::operator - (const Point &p) const {
	assert(dim == p.dim);
	std::vector<double> v = {};
	for (int i = 0; i < dim; i++)
		v.push_back(val[i] - p.val[i]);
	return Point(dim, v);
}

// inner product
double Point::operator * (const Point &p) const {
	assert(dim == p.dim);
	double ret = 0.0;
	for (int i = 0; i < dim; i++)
		ret += val[i] * p.val[i];
	return ret;
}

Point Point::operator * (const double &k) const {
	std::vector<double> v = {};
	for (int i = 0; i < dim; i++) 
		v.push_back(val[i] * k);
	return Point(dim, v);
		
}

// Cross product in two dimension, left X right
double Point::operator ^ (const Point &p) const {
	assert(dim == p.dim);
	assert(dim == 2);
	return val[0] * p.val[1] - val[1] * p.val[0];
}

#pragma endregion

#pragma region Line

Line::Line() {}

Line::~Line() {}

Line::Line(int d, Point a, Point b) {
	dim = d;
	st = a;
	ed = b;
}

int Line::get_dim() {
	return dim;
}

#pragma endregion

#pragma region HalfPlane

HalfPlane::HalfPlane() {}

HalfPlane::~HalfPlane() {}

HalfPlane::HalfPlane(int d, Point &p, double &bx) {
	dim = d;
	a = p;
	b = bx;
}

int HalfPlane::get_dim() {
	return dim;
}

#pragma endregion

#pragma region Polygon

Polygon::Polygon() {}

Polygon::Polygon(int d) {
	dim = d;
}

Polygon::~Polygon() {}

Polygon::Polygon(int d, std::vector<Point> v) {
	vertices.resize(0);
	dim = d;
	for (Point Vi : v)
		vertices.push_back(Vi);
}

int Polygon::get_dim() {
	return dim;
}

#pragma endregion

#pragma region Functions

// Sign function
int sign(const double &x) {
	if (fabs(x) < eps) return 0;
	return x < 0 ? -1 : 1;
}

Point line_intersection(Line &L1, Line &L2) {
	Point u = L1.ed - L1.st, v = L2.ed - L2.st;
	double t1 = ((L2.st - L1.st) ^ v) / (u ^ v);
	return L1.st + u * t1;
}

Line halfplane_boundary(HalfPlane &h) {
	Line ret;
	int dim = h.get_dim();
	if (!sign(h.a.val[0])) {
		std::vector<double> st(dim, 0), ed(dim, 0);
		double val = h.b * (sign(h.a.val[1]) < 0 ? -1 : 1);
		st[0] = 0.0; st[1] = val;
		ed[0] = 1.0; ed[1] = val;
		ret = Line(2, Point(2, st), Point(2, ed));
	}
	else if (!sign(h.a.val[1])) {
		std::vector<double> st(dim, 0), ed(dim, 0);
		double val = h.b * (sign(h.a.val[0]) < 0 ? -1 : 1);
		st[0] = val; st[1] = 0.0;
		ed[0] = val; ed[1] = 1.0;
		ret = Line(2, Point(2, st), Point(2, ed));
	}
	else {
		std::vector<double> st(dim, 0), ed(dim, 0);
		st[0] = 0.0; st[1] = h.b / h.a.val[1];
		ed[0] = 1.0; ed[1] = (h.b - h.a.val[0]) / h.a.val[1];
		ret = Line(2, Point(2, st), Point(2, ed));
	}
	return ret;
}


#pragma endregion

#pragma region LinearProgramming

LinearProgramming::LinearProgramming() {}

LinearProgramming::LinearProgramming(int d) {
	dim = d;
}

LinearProgramming::LinearProgramming(int d, std::vector<HalfPlane> &H_set, Point &obj) {
	dim = d;
	H = H_set;
	O = obj;
}

LinearProgramming::~LinearProgramming() {}

std::vector<Point> LinearProgramming::solve() {

	std::vector<Point> ret = {};

	if (dim == 1) {
		double x_left = -INF, x_right = INF;
		for (int i = 0; i < H.size(); i++) {
			if (sign(H[i].a.val[0]) < 0) 
				x_left = std::max(x_left, -1.0 * H[i].b);
			else if (sign(H[i].a.val[0]) > 0)
				x_right = std::min(x_right, H[i].b);
			else if (sign(H[i].b) < 0) return ret;
		}
		if (sign(x_left - x_right) > 0) return ret;
		else {
			double val_x_left = x_left * O.val[0], val_x_right = x_right * O.val[0];
			std::vector<double> v = {};
			if (sign(val_x_left) >= sign(val_x_right))
				v.push_back(x_left);
			else
				v.push_back(x_right);
			ret.push_back(Point(dim, v));
		}
	}
	else {
		std::vector<double> corner(dim, INF);
		Point v = Point(dim, corner);
		for (int i = 2 * dim; i < H.size(); i++) {
			if (sign(H[i].a * v - H[i].b) <= 0) {
				continue;
			}
			else {
				Line h = halfplane_boundary(H[i]);
				
				std::vector<HalfPlane> new_H = {};
				Point new_O;
				std::vector<double> tmp(1, -1.0);
				Point bound_left = Point(1, tmp);
				tmp[0] = 1.0;
				Point bound_right = Point(1, tmp);

				for (int j = 0; j < i; j++) {
					Line hj = halfplane_boundary(H[j]);
					if (!sign((h.ed - h.st) ^ (hj.ed - hj.st))) {
						if (sign(H[i].a * hj.st - H[i].b) <= 0 || sign(H[j].a * h.st - H[j].b) <= 0) {
							continue;
						}
						else {
							return ret;
						}
					}
					Point inter = line_intersection(h, hj);
					if (sign((inter + (h.ed - h.st) * 0.1) * H[j].a - H[j].b) <= 0) {
						double x = inter.val[0] * -1.0;
						new_H.push_back(HalfPlane(dim - 1, bound_left, x));
					}	
					else
						new_H.push_back(HalfPlane(dim - 1, bound_right, inter.val[0]));	
				}
				if (sign(H[i].a.val[1]) != 0) {
					double x = O.val[0] - O.val[1] * H[i].a.val[0] / H[i].a.val[1];
					std::vector<double> tmp_v = {}; tmp_v.push_back(x);
					new_O = Point(1, tmp_v);
				}
				else {
					new_O = O;
				}

				LinearProgramming new_LP = LinearProgramming(dim - 1, new_H, new_O);
				std::vector<Point> vi = new_LP.solve();
				if (!vi.size()) return ret;
				else {
					std::vector<double> v_tmp(dim, 0);
					v_tmp[0] = vi[0].val[0];
					if (!sign(H[i].a.val[1])) {
						if (sign(H[i].a.val[0]) > 0)
							v_tmp[1] = INF;
						else
							v_tmp[1] = -INF;
					}
					else {
						v_tmp[1] = (H[i].b - H[i].a.val[0] * v_tmp[0]) / H[i].a.val[1];
					}
					v = Point(dim, v_tmp);
				}
			}
		}
		ret.push_back(v);
	}
	return ret;
}

int LinearProgramming::get_dim() {
	return dim;
}

#pragma endregion