#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <math.h>
#include <iostream> 
#include <array>

using namespace std;
namespace py = pybind11;

const double m_ele = 0.51100e-3;
const double m_muon = 0.10566;
const double m_tau = 1.77685;
const double m_pion = 0.13957;

double invertedMatrix(const double m[3][3], double inv[3][3]) {

  double det =
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
    m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
    m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

  if (std::fabs(det)<1e-12) {
    return std::fabs(det);
  }
  
  inv[0][0] =  (m[1][1] * m[2][2] - m[1][2] * m[2][1])/det;
  inv[0][1] = -(m[0][1] * m[2][2] - m[0][2] * m[2][1])/det;
  inv[0][2] =  (m[0][1] * m[1][2] - m[0][2] * m[1][1])/det;

  inv[1][0] = -(m[1][0] * m[2][2] - m[1][2] * m[2][0])/det;
  inv[1][1] =  (m[0][0] * m[2][2] - m[0][2] * m[2][0])/det;
  inv[1][2] = -(m[0][0] * m[1][2] - m[0][2] * m[1][0])/det;
  
  inv[2][0] =  (m[1][0] * m[2][1] - m[1][1] * m[2][0])/det;
  inv[2][1] = -(m[0][0] * m[2][1] - m[0][1] * m[2][0])/det;
  inv[2][2] =  (m[0][0] * m[1][1] - m[0][1] * m[1][0])/det;

  return std::fabs(det);
 
}

vector<double> TauMomGJ(double Evis,
			double cosTheta,
			double mvis) {
  
  double Pvis = sqrt(Evis*Evis-mvis*mvis);

  vector<double> solutions;
  
  double A = 4 * ( Evis*Evis - Pvis*Pvis*cosTheta*cosTheta );
  double B = 4 * Pvis * ( m_tau*m_tau + mvis*mvis) * cosTheta;
  double C = 4 * Evis*Evis*m_tau*m_tau - (m_tau*m_tau + mvis*mvis)*(m_tau*m_tau + mvis*mvis);
  double disc = B*B - 4*A*C;

  if (disc<0) {
    solutions.push_back(-1.0);
    return solutions;
  }
  
  double result = (B - sqrt(disc)) / (2 * A);
  solutions.push_back(result);

  result = (B + sqrt(disc)) / (2 * A);
  solutions.push_back(result);

  return solutions;
  
}

double Dot(const double v1[3],
	   const double v2[3]) {
  
  double dot = 0;
  for (int i=0; i<3; ++i) {
    dot += v1[i]*v2[i];
  }
  return dot;
  
}

void Cross(const double v1[3],
	   const double v2[3],
	   double out[3]) {

  out[0] = v1[1]*v2[2] - v1[2]*v2[1];
  out[1] = v1[2]*v2[0] - v1[0]*v2[2];
  out[2] = v1[0]*v2[1] - v1[1]*v2[0];

}

void LinearSum(const double a,
	       const double b,
	       const double va[3],
	       const double vb[3],
	       double out[3]) {
  
  for (int i=0; i<3; ++i) {
    out[i] = a*va[i] + b*vb[i];    
  }
  
}

void LinearSum3(const double a,
		const double b,
		const double c,
		const double va[3],
		const double vb[3],
		const double vc[3],
		double out[3]) {
  
  for (int i=0; i<3; ++i) {
    out[i] = a*va[i] + b*vb[i] + c*vc[i];    
  }
  
}

double Convolution(const double m[3][3], const double a[3]) {
  
  double x = 0;
  for (int i=0; i<3; ++i) {
    for (int j=0; j<3; ++j) {
      x += a[i]*m[i][j]*a[j];
    }
  }
  return x;
  
}

double ThetaGJ(double Pvis, double Ptau, double mvis) {

  double Evis = sqrt(Pvis*Pvis+mvis*mvis);
  double Etau = sqrt(Ptau*Ptau+m_tau*m_tau);
  double numerator = 2*Evis*Etau - m_tau*m_tau - mvis*mvis;
  double denominator = 2*Ptau*Pvis;
  double cosTheta = max(-1.,min(1.,numerator/denominator));
  double theta = acos(cosTheta);
  return theta;

}

double ThetaGJmax(double Pvis, double mvis) {

  double sinTheta = 0.5*(m_tau*m_tau-mvis*mvis)/(m_tau*Pvis);
  if (sinTheta<-0.99) sinTheta = -0.99;
  if (sinTheta>0.99) sinTheta = 0.99;
  double theta = asin(sinTheta);
  return theta;
  
}
