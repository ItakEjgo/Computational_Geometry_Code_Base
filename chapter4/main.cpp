#include <iostream>
#include <vector>
#include "geo_header.h"
#include <fstream>
#include <algorithm>

using namespace std;

vector<HalfPlane> H;


int main()
{
	ifstream fin("input1.txt");

	int d, n;
	fin >> d >> n;

	for (int i = 0; i < d; i++) {
		vector<double> a(d, 0);
		a[i] = 1.0;
		double b = INF;
		Point p = Point(d, a);
		HalfPlane hp = HalfPlane(d, p, b);
		H.push_back(hp);
	}

	for (int i = 0; i < d; i++) {
		vector<double> a(d, 0);
		a[i] = -1.0;
		double b = INF;
		Point p = Point(d, a);
		HalfPlane hp = HalfPlane(d, p, b);
		H.push_back(hp);
	}

	for (int i = 0; i < n; i++) {
		vector<double> a(d, 0);
		double b;
		for (int j = 0; j < d; j++) 
			fin >> a[j];
		fin >> b;
		Point p = Point(d, a);
		HalfPlane hp = HalfPlane(d, p, b);
		H.push_back(hp);
	}
	
	vector<double> a(d, 0);
	for (int i = 0; i < d; i++) {
		fin >> a[i];
	}
	Point O = Point(d, a);
	
	// random_shuffle(H.begin() + 2 * d, H.end());

	LinearProgramming LP = LinearProgramming(d, H, O);
	vector<Point> res = LP.solve();
	if (!res.size()) {
		cout << "Infeasible" << endl;
	}
	for (int i = 0; i < res.size(); i++){
		for (int j = 0; j < d; j++) {
			cout << res[i].val[j] << " ";
		}
		cout << res[i] * O << endl;
	}

	return 0; 
}

