#ifndef _COMMON_H_
#define _COMMON_H_

#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <array>

using std::vector;
using std::map;
using std::pair;
using std::array;

#define ERROR(msg) do {fprintf(stderr, "[E::%s] %s\n", __func__, msg); ::abort();} while(0)

//class Coord
//{
//public:
//  explicit Coord():
//    m_x(0), m_y(0), m_z(0)
//  {}
//  explicit Coord(double x, double y, double z):
//    m_x(x), m_y(y), m_z(z)
//  {}
//  double m_x;
//  double m_y;
//  double m_z;
//
//  static double distance(Coord a, Coord b)
//  {
//    return (a.m_x-b.m_x)*(a.m_x - b.m_x)
//      + (a.m_x-b.m_x)*(a.m_x - b.m_x)
//      + (a.m_x-b.m_x)*(a.m_x - b.m_x);
//  }
//};

double distance(const vector<double>& vec1, const vector<double>& vec2);

vector<double> translate(const vector<double>& vec, const vector<double>& dvec);
double dot(const vector<double>& v1, const vector<double>&v2);
vector<double> div(const vector<double>& vec, double factor);
vector<double> rotate(const vector<double>& vec, const vector<double>& axis, double theta);
vector<double> get_normal(const vector<double>& vec1, const vector<double>& vec2);
double length(const vector<double>& vec);
double angle(const vector<double>& vec1, const vector<double>& vec2);

double linear1(double dxx0, double dx1x0, double y0, double y1);
double bilinear(double y0, double y1, double y2, double y3, double angx, double angy);
double bilinear_gen(double y0, double y1, double y2, double y3, double angx1, double angx2, double angy, int label=1);

double lagrange_interp(const std::vector<std::pair<double, double>>& points, double x);

// return a-b
vector<double> subtraction(const vector<double>& a, const vector<double>& b);
vector<double> get_bisect_unit(const vector<double>& a, const vector<double>& b);
vector<double> get_normal_unit(const vector<double>& a, const vector<double>& b);
void get_unit(vector<double>& vec);

/**
 *   Return the weights of the two configurations at one bisector direction.
 *   Linear interpolation is used.
 *   NOTE: The two vectors of the two configurations are calculated here
 *       but in future, they should be retrieve from a pre-calculated
 *       constant data list or dictionary of fixed directions.
 */
vector<int> weights_in_subsection(const vector<double>& bisvec, double& wghx, double& wghy, double cutoff=0.9999);

void weights_for_normal_general(const vector<double>& normal_vec, const vector<vector<double>>& config_vecs, double& w1, double& w2, int& ndx1, int& ndx2, double cutoff=0.0000001);

pair<double, array<int, 3>> get_neighors_for_normal(const vector<double>& normal_vec, const vector<vector<double>>& config_vecs, double cutoff=0.0000001);

#endif
