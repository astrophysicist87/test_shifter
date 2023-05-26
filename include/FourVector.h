#ifndef SHIFT_LIB_FOURVECTOR_H
#define SHIFT_LIB_FOURVECTOR_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>

using namespace std;

namespace shift_lib
{


	// Powers of small integers - for balance speed/code clarity.
	inline double pow2(const double& x) {return x*x;}
	inline double pow3(const double& x) {return x*x*x;}
	inline double pow4(const double& x) {return x*x*x*x;}
	inline double pow5(const double& x) {return x*x*x*x*x;}
	inline double pow6(const double& x) {return x*x*x*x*x*x;}
	inline double pow7(const double& x) {return x*x*x*x*x*x*x;}
	inline double pow8(const double& x) {return x*x*x*x*x*x*x*x;}

	// Avoid problem with negative square root argument (from roundoff).
	inline double sqrtpos(const double& x) {return sqrt( max( 0., x));}


	class Vec4 {

	public:

	  // Constructors.
	  Vec4(double xIn = 0., double yIn = 0., double zIn = 0., double tIn = 0.)
		: xx(xIn), yy(yIn), zz(zIn), tt(tIn) { }
	  Vec4(const Vec4& v) : xx(v.xx), yy(v.yy), zz(v.zz), tt(v.tt) { }
	  Vec4& operator=(const Vec4& v) { if (this != &v) { xx = v.xx; yy = v.yy;
		zz = v.zz; tt = v.tt; } return *this; }
	  Vec4& operator=(double value) { xx = value; yy = value; zz = value;
		tt = value; return *this; }

	  // Member functions for input.
	  void reset() {xx = 0.; yy = 0.; zz = 0.; tt = 0.;}
	  void p(double xIn, double yIn, double zIn, double tIn)
		{xx = xIn; yy = yIn; zz = zIn; tt = tIn;}
	  void p(Vec4 pIn) {xx = pIn.xx; yy = pIn.yy; zz = pIn.zz; tt = pIn.tt;}
	  void x(double xIn) {xx = xIn;}
	  void y(double yIn) {yy = yIn;}
	  void z(double zIn) {zz = zIn;}
	  void t(double tIn) {tt = tIn;}

	  // Member functions for output.
	  double x() const {return xx;}
	  double y() const {return yy;}
	  double z() const {return zz;}
	  double t() const {return tt;}
	  double& operator[](int i) {
		if      (i == 1) return xx;
		else if (i == 2) return yy;
		else if (i == 3) return zz;
		else             return tt;
	  }
	  double mCalc() const {double temp = tt*tt - xx*xx - yy*yy - zz*zz;
		return (temp >= 0.) ? sqrt(temp) : -sqrt(-temp);}
	  double m2Calc() const {return tt*tt - xx*xx - yy*yy - zz*zz;}
	  double pT() const {return sqrt(xx*xx + yy*yy);}
	  double pT2() const {return xx*xx + yy*yy;}
	  double pAbs() const {return sqrt(xx*xx + yy*yy + zz*zz);}
	  double pAbs2() const {return xx*xx + yy*yy + zz*zz;}
	  double eT() const {double temp = xx*xx + yy*yy;
		return tt * sqrt( temp / (temp + zz*zz) );}
	  double eT2() const {double temp = xx*xx + yy*yy;
		return tt*tt * temp / (temp + zz*zz);}
	  double theta() const {return atan2(sqrt(xx*xx + yy*yy), zz);}
	  double phi() const {return atan2(yy,xx);}
	  double thetaXZ() const {return atan2(xx,zz);}
	  double pPos() const {return tt + zz;}
	  double pNeg() const {return tt - zz;}
	  double rap() const {return 0.5 * log( (tt + zz) / (tt - zz) );}
	  double eta() const {double xyz = sqrt(xx*xx + yy*yy + zz*zz);
		return 0.5 * log( (xyz + zz) / (xyz - zz) );}

	  // Member functions that perform operations.
	  void rescale3(double fac) {xx *= fac; yy *= fac; zz *= fac;}
	  void rescale4(double fac) {xx *= fac; yy *= fac; zz *= fac; tt *= fac;}
	  void flip3() {xx = -xx; yy = -yy; zz = -zz;}
	  void flip4() {xx = -xx; yy = -yy; zz = -zz; tt = -tt;}
	  void rot(double thetaIn, double phiIn);
	  void rotaxis(double phiIn, double nx, double ny, double nz);
	  void rotaxis(double phiIn, const Vec4& n);
	  void bst(double betaX, double betaY, double betaZ);
	  void bst(double betaX, double betaY, double betaZ, double gamma);
	  void bst(const Vec4& pIn);
	  void bst(const Vec4& pIn, double mIn);
	  void bstback(const Vec4& pIn);
	  void bstback(const Vec4& pIn, double mIn);

	  // Operator overloading with member functions
	  inline Vec4 operator-() const {Vec4 tmp; tmp.xx = -xx; tmp.yy = -yy;
		tmp.zz = -zz; tmp.tt = -tt; return tmp;}
	  inline Vec4& operator+=(const Vec4& v) {xx += v.xx; yy += v.yy; zz += v.zz;
		tt += v.tt; return *this;}
	  inline Vec4& operator-=(const Vec4& v) {xx -= v.xx; yy -= v.yy; zz -= v.zz;
		tt -= v.tt; return *this;}
	  inline Vec4& operator*=(double f) {xx *= f; yy *= f; zz *= f;
		tt *= f; return *this;}
	  inline Vec4& operator/=(double f) {xx /= f; yy /= f; zz /= f;
		tt /= f; return *this;}
	  inline Vec4 operator+(const Vec4& v) const {
		Vec4 tmp = *this; return tmp += v;}
	  inline Vec4 operator-(const Vec4& v) const {
		Vec4 tmp = *this; return tmp -= v;}
	  inline Vec4 operator*(double f) const {
		Vec4 tmp = *this; return tmp *= f;}
	  inline Vec4 operator/(double f) const {
		Vec4 tmp = *this; return tmp /= f;}
	  inline double operator*(const Vec4& v) const {
		return tt*v.tt - xx*v.xx - yy*v.yy - zz*v.zz;}

	  // Operator overloading with friends
	  friend Vec4 operator*(double f, const Vec4& v1);

	  // Print a four-vector.
	  friend ostream& operator<<(ostream&, const Vec4& v) ;

	  // Invariant mass of a pair and its square.
	  friend double m(const Vec4& v1, const Vec4& v2);
	  friend double m2(const Vec4& v1, const Vec4& v2);

	  // Scalar and cross product of 3-vector parts.
	  friend double dot3(const Vec4& v1, const Vec4& v2);
	  friend Vec4 cross3(const Vec4& v1, const Vec4& v2);

	  // Cross-product of three 4-vectors ( p_i = epsilon_{iabc} p_a p_b p_c).
	  friend Vec4 cross4(const Vec4& a, const Vec4& b, const Vec4& c);

	  // theta is polar angle between v1 and v2.
	  friend double theta(const Vec4& v1, const Vec4& v2);
	  friend double costheta(const Vec4& v1, const Vec4& v2);

	  // phi is azimuthal angle between v1 and v2 around z axis.
	  friend double phi(const Vec4& v1, const Vec4& v2);
	  friend double cosphi(const Vec4& v1, const Vec4& v2);

	  // phi is azimuthal angle between v1 and v2 around n axis.
	  friend double phi(const Vec4& v1, const Vec4& v2, const Vec4& n);
	  friend double cosphi(const Vec4& v1, const Vec4& v2, const Vec4& n);

	  // R is distance in cylindrical (y/eta, phi) coordinates.
	  friend double RRapPhi(const Vec4& v1, const Vec4& v2);
	  friend double REtaPhi(const Vec4& v1, const Vec4& v2);

	  // Shift four-momenta within pair from old to new masses.
	  friend bool pShift( Vec4& p1Move, Vec4& p2Move, double m1New, double m2New);

	  // Create two vectors that are perpendicular to both input vectors.
	  friend pair<Vec4,Vec4> getTwoPerpendicular(const Vec4& v1, const Vec4& v2);

	private:

	  // Constants: could only be changed in the code itself.
	  static const double TINY;

	  // The four-vector data members.
	  double xx, yy, zz, tt;

	};

	//--------------------------------------------------------------------------

	// Namespace function declarations; friends of Vec4 class.

	// Implementation of operator overloading with friends.
	inline Vec4 operator*(double f, const Vec4& v1)
	  {Vec4 v = v1; return v *= f;}

	// Invariant mass of a pair and its square.
	double m(const Vec4& v1, const Vec4& v2);
	double m2(const Vec4& v1, const Vec4& v2);

	// Scalar and cross product of 3-vector parts.
	double dot3(const Vec4& v1, const Vec4& v2);
	Vec4 cross3(const Vec4& v1, const Vec4& v2);

	// Cross-product of three 4-vectors ( p_i = epsilon_{iabc} p_a p_b p_c).
	Vec4 cross4(const Vec4& a, const Vec4& b, const Vec4& c);

	// theta is polar angle between v1 and v2.
	double theta(const Vec4& v1, const Vec4& v2);
	double costheta(const Vec4& v1, const Vec4& v2);

	// phi is azimuthal angle between v1 and v2 around z axis.
	double phi(const Vec4& v1, const Vec4& v2);
	double cosphi(const Vec4& v1, const Vec4& v2);

	// phi is azimuthal angle between v1 and v2 around n axis.
	double phi(const Vec4& v1, const Vec4& v2, const Vec4& n);
	double cosphi(const Vec4& v1, const Vec4& v2, const Vec4& n);

	// R is distance in cylindrical (y/eta, phi) coordinates.
	double RRapPhi(const Vec4& v1, const Vec4& v2);
	double REtaPhi(const Vec4& v1, const Vec4& v2);

	// Print a four-vector.
	ostream& operator<<(ostream&, const Vec4& v) ;

	// Shift four-momenta within pair from old to new masses.
	bool pShift( Vec4& p1Move, Vec4& p2Move, double m1New, double m2New);

	// Create two vectors that are perpendicular to both input vectors.
	pair<Vec4,Vec4> getTwoPerpendicular(const Vec4& v1, const Vec4& v2);

}

#endif
