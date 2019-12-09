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
    double operator ^ (const Point &p) const { return x * p.y - y * p.x; }
    double operator * (const Point &p) const { return x * p.x + y * p.y; }
    bool operator == (const Point &p) const { return !sign(x - p.x) && !sign(y - p.y);}
};

typedef vector<Point> VP;
typedef pair<int, int> pii;
typedef pair<double, double> pdd;

Point origin;
VP P;

//check whether r lies left of line pq
bool lies_left(Point p, Point q, Point r) {
    return sign((q - p) ^ (r - p)) > 0 ? 1 : 0;
}

bool cmp(const Point &a, const Point &b) {
    //return sign(angle(a) - angle(b)) < 0;
    return sign((a - origin) ^ (b - origin)) < 0;
}

VP slow_convex_hull(VP P){
    VP ret = {};
    for (int i = 0; i < P.size(); i++){
        for (int j = 0; j < P.size(); j++){ //  enumerate all pairs
            if (i == j) continue;
            bool valid = 1;
            for (int r = 0; r < P.size(); r++){ //  check whether all points in the same side
                if (r == i || r == j) continue;
                if (lies_left(P[i], P[j], P[r])) {
                    valid = 0;
                    break;
                }
            }
            if (valid) ret.push_back(P[i]);
        }
    }
    origin = ret[0];
    sort(ret.begin(), ret.end(), cmp);  //  sort points in counter-clockwise order
    VP res = {}; 
    res.push_back(ret[0]);
    for (int i = 1; i < ret.size(); i++){
       if (!(ret[i] == ret[i - 1])) res.push_back(ret[i]);  // remove duplicates 
    }
    return res;
}

int main(){
    int n; cin >> n;
    for (int i = 0; i < n; ++i){
        double x, y; cin >> x >> y;
        P.push_back(Point(x, y));
    }
    VP res = slow_convex_hull(P);
    for (int i = 0; i < res.size(); i++){
        cout << res[i].x << " " << res[i].y << endl;
    }
    return 0;
}
