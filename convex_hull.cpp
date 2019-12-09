#include "bits/stdc++.h"

using namespace std;

const double eps = 1e-9;
int sign(const double x){
    if (fabs(x) < eps) return 0;
    return x < 0 ? -1 : 1;
}

struct Point{
    double x, y;
    Point(){}
    Point(double x, double y):x(x), y(y){}
    Point operator + (const Point &p) const { return Point(x + p.x, y + p.y); }
    Point operator - (const Point &p) const { return Point(x - p.x, y - p.y); }
    double operator ^ (const Point &p) const { return x * p.y - y * p.x; }  //  cross product
    double operator * (const Point &p) const { return x * p.x + y * p.y; }  //  inner product
    bool operator < (const Point &p) const{ return x < p.x || x == p.x && y < p.y; }
};

typedef vector<Point> VP;
typedef pair<int, int> pii;
typedef pair<double, double> pdd;

VP P;

//check whether r lies left or on the line pq
bool lies_left(Point p, Point q, Point r) {
    return sign((q - p) ^ (r - p)) > 0 ? 1 : 0;
}

bool cmp(const pii &a, const pii &b) {
    return sign((P[a.second] - P[a.first]) ^ (P[b.second] - P[a.first])) <= 0;
}

VP convex_hull(VP P){
    sort(P.begin(), P.end());  //   sort all points
    VP ret = {};
    int sz = 0;
    for (int i = 0; i < P.size(); i++) { // upper hull
        while (sz > 1 && !lies_left(ret[sz - 2], P[i], ret[sz - 1])) {
            ret.pop_back();
            sz--;
        }
        ret.push_back(P[i]);
        sz++;
    }
    int k = sz;
    for (int i = P.size() - 2; i >= 0; i--){ // lower hull
        while (sz > k && !lies_left(ret[sz - 2], P[i], ret[sz - 1])) {
            ret.pop_back();
            sz--;
        }
        ret.push_back(P[i]);
        sz++;
    }
    ret.pop_back();
    return ret;
}

int main(){
    int n; cin >> n;
    for (int i = 0; i < n; ++i){
        double x, y; cin >> x >> y;
        P.push_back(Point(x, y));
    }
    VP res = convex_hull(P);
    for (int i = 0; i < res.size(); i++){
        cout << res[i].x << " " << res[i].y << endl;
    }
    return 0;
}
