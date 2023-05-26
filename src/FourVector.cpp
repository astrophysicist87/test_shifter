#include "../include/FourVector.h"

namespace shift_lib
{


	//==========================================================================

	// Vec4 class.
	// This class implements four-vectors, in energy-momentum space.
	// (But could also be used to hold space-time four-vectors.)

	//--------------------------------------------------------------------------

	// Constants: could be changed here if desired, but normally should not.
	// These are of technical nature, as described for each.

	// Small number to avoid division by zero.
	const double Vec4::TINY = 1e-20;

	//--------------------------------------------------------------------------

	// Rotation (simple).

	void Vec4::rot(double thetaIn, double phiIn) {

	  double cthe = cos(thetaIn);
	  double sthe = sin(thetaIn);
	  double cphi = cos(phiIn);
	  double sphi = sin(phiIn);
	  double tmpx =  cthe * cphi * xx -    sphi * yy + sthe * cphi * zz;
	  double tmpy =  cthe * sphi * xx +    cphi * yy + sthe * sphi * zz;
	  double tmpz = -sthe *        xx +                cthe *        zz;
	  xx          = tmpx;
	  yy          = tmpy;
	  zz          = tmpz;

	}

	//--------------------------------------------------------------------------

	// Azimuthal rotation phi around an arbitrary axis (nz, ny, nz).

	void Vec4::rotaxis(double phiIn, double nx, double ny, double nz) {

	  double norm = 1./sqrt(nx*nx + ny*ny + nz*nz);
	  nx         *= norm;
	  ny         *= norm;
	  nz         *= norm;
	  double cphi = cos(phiIn);
	  double sphi = sin(phiIn);
	  double comb = (nx * xx + ny * yy + nz * zz) * (1. - cphi);
	  double tmpx = cphi * xx + comb * nx + sphi * (ny * zz - nz * yy);
	  double tmpy = cphi * yy + comb * ny + sphi * (nz * xx - nx * zz);
	  double tmpz = cphi * zz + comb * nz + sphi * (nx * yy - ny * xx);
	  xx          = tmpx;
	  yy          = tmpy;
	  zz          = tmpz;

	}

	//--------------------------------------------------------------------------

	// Azimuthal rotation phi around an arbitrary (3-vector component of) axis.

	void Vec4::rotaxis(double phiIn, const Vec4& n) {

	  double nx   = n.xx;
	  double ny   = n.yy;
	  double nz   = n.zz;
	  double norm = 1./sqrt(nx*nx + ny*ny + nz*nz);
	  nx         *= norm;
	  ny          *=norm;
	  nz          *=norm;
	  double cphi = cos(phiIn);
	  double sphi = sin(phiIn);
	  double comb = (nx * xx + ny * yy + nz * zz) * (1. - cphi);
	  double tmpx = cphi * xx + comb * nx + sphi * (ny * zz - nz * yy);
	  double tmpy = cphi * yy + comb * ny + sphi * (nz * xx - nx * zz);
	  double tmpz = cphi * zz + comb * nz + sphi * (nx * yy - ny * xx);
	  xx          = tmpx;
	  yy          = tmpy;
	  zz          = tmpz;

	}

	//--------------------------------------------------------------------------

	// Boost (simple).

	void Vec4::bst(double betaX, double betaY, double betaZ) {

	  double beta2 = betaX*betaX + betaY*betaY + betaZ*betaZ;
	  if (beta2 >= 1.) return;
	  double gamma = 1. / sqrt(1. - beta2);
	  double prod1 = betaX * xx + betaY * yy + betaZ * zz;
	  double prod2 = gamma * (gamma * prod1 / (1. + gamma) + tt);
	  xx += prod2 * betaX;
	  yy += prod2 * betaY;
	  zz += prod2 * betaZ;
	  tt = gamma * (tt + prod1);

	}

	//--------------------------------------------------------------------------

	// Boost (simple, given gamma).

	void Vec4::bst(double betaX, double betaY, double betaZ, double gamma) {

	  double prod1 = betaX * xx + betaY * yy + betaZ * zz;
	  double prod2 = gamma * (gamma * prod1 / (1. + gamma) + tt);
	  xx += prod2 * betaX;
	  yy += prod2 * betaY;
	  zz += prod2 * betaZ;
	  tt = gamma * (tt + prod1);

	}

	//--------------------------------------------------------------------------

	// Boost given by a Vec4 p.

	void Vec4::bst(const Vec4& pIn) {

	  if (abs(pIn.tt) < Vec4::TINY) return;
	  double betaX = pIn.xx / pIn.tt;
	  double betaY = pIn.yy / pIn.tt;
	  double betaZ = pIn.zz / pIn.tt;
	  double beta2 = betaX*betaX + betaY*betaY + betaZ*betaZ;
	  if (beta2 >= 1.) return;
	  double gamma = 1. / sqrt(1. - beta2);
	  double prod1 = betaX * xx + betaY * yy + betaZ * zz;
	  double prod2 = gamma * (gamma * prod1 / (1. + gamma) + tt);
	  xx          += prod2 * betaX;
	  yy          += prod2 * betaY;
	  zz          += prod2 * betaZ;
	  tt           = gamma * (tt + prod1);

	}

	//--------------------------------------------------------------------------

	// Boost given by a Vec4 p and double m.

	void Vec4::bst(const Vec4& pIn, double mIn) {

	  if (abs(pIn.tt) < Vec4::TINY) return;
	  double betaX = pIn.xx / pIn.tt;
	  double betaY = pIn.yy / pIn.tt;
	  double betaZ = pIn.zz / pIn.tt;
	  double gamma = pIn.tt / mIn;
	  double prod1 = betaX * xx + betaY * yy + betaZ * zz;
	  double prod2 = gamma * (gamma * prod1 / (1. + gamma) + tt);
	  xx          += prod2 * betaX;
	  yy          += prod2 * betaY;
	  zz          += prod2 * betaZ;
	  tt           = gamma * (tt + prod1);

	}

	//--------------------------------------------------------------------------

	// Boost given by a Vec4 p; boost in opposite direction.

	void Vec4::bstback(const Vec4& pIn) {

	  if (abs(pIn.tt) < Vec4::TINY) return;
	  double betaX = -pIn.xx / pIn.tt;
	  double betaY = -pIn.yy / pIn.tt;
	  double betaZ = -pIn.zz / pIn.tt;
	  double beta2 = betaX*betaX + betaY*betaY + betaZ*betaZ;
	  if (beta2 >= 1.) return;
	  double gamma = 1. / sqrt(1. - beta2);
	  double prod1 = betaX * xx + betaY * yy + betaZ * zz;
	  double prod2 = gamma * (gamma * prod1 / (1. + gamma) + tt);
	  xx          += prod2 * betaX;
	  yy          += prod2 * betaY;
	  zz          += prod2 * betaZ;
	  tt           = gamma * (tt + prod1);

	}

	//--------------------------------------------------------------------------

	// Boost given by a Vec4 p and double m; boost in opposite direction.

	void Vec4::bstback(const Vec4& pIn, double mIn) {

	  if (abs(pIn.tt) < Vec4::TINY) return;
	  double betaX = -pIn.xx / pIn.tt;
	  double betaY = -pIn.yy / pIn.tt;
	  double betaZ = -pIn.zz / pIn.tt;
	  double gamma = pIn.tt / mIn;
	  double prod1 = betaX * xx + betaY * yy + betaZ * zz;
	  double prod2 = gamma * (gamma * prod1 / (1. + gamma) + tt);
	  xx          += prod2 * betaX;
	  yy          += prod2 * betaY;
	  zz          += prod2 * betaZ;
	  tt           = gamma * (tt + prod1);

	}


	//--------------------------------------------------------------------------

	// Print a four-vector: also operator overloading with friend.

	ostream& operator<<(ostream& os, const Vec4& v) {
	  os << fixed << setprecision(3) << " " << setw(9) << v.xx << " "
		 << setw(9) << v.yy << " " << setw(9) << v.zz << " " << setw(9)
		 << v.tt << " (" << setw(9) << v.mCalc() << ")\n";
	  return os;
	}

	//--------------------------------------------------------------------------

	// The invariant mass of two four-vectors.

	double m(const Vec4& v1, const Vec4& v2) {
	  double m2 = pow2(v1.tt + v2.tt) - pow2(v1.xx + v2.xx)
		 - pow2(v1.yy + v2.yy) - pow2(v1.zz + v2.zz);
	  return (m2 > 0.) ? sqrt(m2) : 0.;
	}

	//--------------------------------------------------------------------------

	// The squared invariant mass of two four-vectors.

	double m2(const Vec4& v1, const Vec4& v2) {
	  double m2 = pow2(v1.tt + v2.tt) - pow2(v1.xx + v2.xx)
		 - pow2(v1.yy + v2.yy) - pow2(v1.zz + v2.zz);
	  return m2;
	}

	//--------------------------------------------------------------------------

	// The scalar product of two three-vectors.

	double dot3(const Vec4& v1, const Vec4& v2) {
	  return v1.xx*v2.xx + v1.yy*v2.yy + v1.zz*v2.zz;
	}

	//--------------------------------------------------------------------------

	// The cross product of two three-vectors.

	Vec4 cross3(const Vec4& v1, const Vec4& v2) {
	  Vec4 v;
	  v.xx = v1.yy * v2.zz - v1.zz * v2.yy;
	  v.yy = v1.zz * v2.xx - v1.xx * v2.zz;
	  v.zz = v1.xx * v2.yy - v1.yy * v2.xx; return v;
	}


	//--------------------------------------------------------------------------

	// Cross-product of three 4-vectors ( p_i = epsilon_{iabc} p_a p_b p_c)

	Vec4 cross4(const Vec4& a, const Vec4& b, const Vec4& c) {
	  Vec4 v(0.,0.,0.,0.);
	  v.tt =   a.xx*b.yy*c.zz + a.yy*b.zz*c.xx + a.zz*b.xx*c.yy
		     - a.xx*b.zz*c.yy - a.zz*b.yy*c.xx - a.yy*b.xx*c.zz;
	  v.xx = -(- a.tt*b.yy*c.zz - a.yy*b.zz*c.tt - a.zz*b.tt*c.yy
		       + a.tt*b.zz*c.yy + a.zz*b.yy*c.tt + a.yy*b.tt*c.zz);
	  v.yy = -(- a.xx*b.tt*c.zz - a.tt*b.zz*c.xx - a.zz*b.xx*c.tt
		       + a.xx*b.zz*c.tt + a.zz*b.tt*c.xx + a.tt*b.xx*c.zz);
	  v.zz = -(- a.xx*b.yy*c.tt - a.yy*b.tt*c.xx - a.tt*b.xx*c.yy
		       + a.xx*b.tt*c.yy + a.tt*b.yy*c.xx + a.yy*b.xx*c.tt);
	  return v;
	}

	//--------------------------------------------------------------------------

	// Opening angle between two three-vectors.

	double theta(const Vec4& v1, const Vec4& v2) {
	  double cthe = (v1.xx * v2.xx + v1.yy * v2.yy + v1.zz * v2.zz)
		/ sqrt( (v1.xx*v1.xx + v1.yy*v1.yy + v1.zz*v1.zz)
		* (v2.xx*v2.xx + v2.yy*v2.yy + v2.zz*v2.zz) );
	  cthe = max(-1., min(1., cthe));
	  return acos(cthe);
	}

	//--------------------------------------------------------------------------

	// Cosine of the opening angle between two three-vectors.

	double costheta(const Vec4& v1, const Vec4& v2) {
	  double cthe = (v1.xx * v2.xx + v1.yy * v2.yy + v1.zz * v2.zz)
		/ sqrt( (v1.xx*v1.xx + v1.yy*v1.yy + v1.zz*v1.zz)
		* (v2.xx*v2.xx + v2.yy*v2.yy + v2.zz*v2.zz) );
	  cthe = max(-1., min(1., cthe));
	  return cthe;
	}

	//--------------------------------------------------------------------------

	// Azimuthal angle between two three-vectors.

	double phi(const Vec4& v1, const Vec4& v2) {
	  double cphi = (v1.xx * v2.xx + v1.yy * v2.yy) / sqrt( max( Vec4::TINY,
		(v1.xx*v1.xx + v1.yy*v1.yy) * (v2.xx*v2.xx + v2.yy*v2.yy) ));
	  cphi = max(-1., min(1., cphi));
	  return acos(cphi);
	}

	//--------------------------------------------------------------------------

	// Cosine of the azimuthal angle between two three-vectors.

	double cosphi(const Vec4& v1, const Vec4& v2) {
	  double cphi = (v1.xx * v2.xx + v1.yy * v2.yy) / sqrt( max( Vec4::TINY,
		(v1.xx*v1.xx + v1.yy*v1.yy) * (v2.xx*v2.xx + v2.yy*v2.yy) ));
	  cphi = max(-1., min(1., cphi));
	  return cphi;
	}

	//--------------------------------------------------------------------------

	// Azimuthal angle between two three-vectors around a third.

	double phi(const Vec4& v1, const Vec4& v2, const Vec4& n) {
	  double nx = n.xx; double ny = n.yy; double nz = n.zz;
	  double norm = 1. / sqrt(nx*nx + ny*ny + nz*nz);
	  nx *= norm; ny *=norm; nz *=norm;
	  double v1s = v1.xx * v1.xx + v1.yy * v1.yy + v1.zz * v1.zz;
	  double v2s = v2.xx * v2.xx + v2.yy * v2.yy + v2.zz * v2.zz;
	  double v1v2 = v1.xx * v2.xx + v1.yy * v2.yy + v1.zz * v2.zz;
	  double v1n = v1.xx * nx + v1.yy * ny + v1.zz * nz;
	  double v2n = v2.xx * nx + v2.yy * ny + v2.zz * nz;
	  double cphi = (v1v2 - v1n * v2n) / sqrt( max( Vec4::TINY,
		(v1s - v1n*v1n) * (v2s - v2n*v2n) ));
	  cphi = max(-1., min(1., cphi));
	  return acos(cphi);
	}

	//--------------------------------------------------------------------------

	// Cosine of the azimuthal angle between two three-vectors around a third.

	double cosphi(const Vec4& v1, const Vec4& v2, const Vec4& n) {
	  double nx = n.xx; double ny = n.yy; double nz = n.zz;
	  double norm = 1. / sqrt(nx*nx + ny*ny + nz*nz);
	  nx *= norm; ny *=norm; nz *=norm;
	  double v1s = v1.xx * v1.xx + v1.yy * v1.yy + v1.zz * v1.zz;
	  double v2s = v2.xx * v2.xx + v2.yy * v2.yy + v2.zz * v2.zz;
	  double v1v2 = v1.xx * v2.xx + v1.yy * v2.yy + v1.zz * v2.zz;
	  double v1n = v1.xx * nx + v1.yy * ny + v1.zz * nz;
	  double v2n = v2.xx * nx + v2.yy * ny + v2.zz * nz;
	  double cphi = (v1v2 - v1n * v2n) / sqrt( max( Vec4::TINY,
		(v1s - v1n*v1n) * (v2s - v2n*v2n) ));
	  cphi = max(-1., min(1., cphi));
	  return cphi;
	}

	//--------------------------------------------------------------------------

	// Distance in cylindrical (y, phi) coordinates.

	double RRapPhi(const Vec4& v1, const Vec4& v2) {
	  double dRap = abs(v1.rap() - v2.rap());
	  double dPhi = abs(v1.phi() - v2.phi());
	  if (dPhi > M_PI) dPhi = 2. * M_PI - dPhi;
	  return sqrt(dRap*dRap + dPhi*dPhi);
	}

	//--------------------------------------------------------------------------

	// Distance in cylindrical (eta, phi) coordinates.

	double REtaPhi(const Vec4& v1, const Vec4& v2) {
	  double dEta = abs(v1.eta() - v2.eta());
	  double dPhi = abs(v1.phi() - v2.phi());
	  if (dPhi > M_PI) dPhi = 2. * M_PI - dPhi;
	  return sqrt(dEta*dEta + dPhi*dPhi);
	}

	//--------------------------------------------------------------------------

	// Shift four-momenta within pair from old to new masses.
	// Note that p1Move and p2Move change values during operation.

	bool pShift( Vec4& p1Move, Vec4& p2Move, double m1New, double m2New) {

	  // Standard kinematics variables.
	  double sH  = (p1Move + p2Move).m2Calc();
	  double r1  = p1Move.m2Calc() / sH;
	  double r2  = p2Move.m2Calc() / sH;
	  double r3  = m1New * m1New / sH;
	  double r4  = m2New * m2New / sH;
	  double l12 = sqrtpos(pow2(1. - r1 - r2) - 4. * r1 * r2);
	  double l34 = sqrtpos(pow2(1. - r3 - r4) - 4. * r3 * r4);

	  // Check that shift operation possible.
	  if (sH <= pow2(m1New + m2New) || l12 < Vec4::TINY || l34 < Vec4::TINY)
		return false;

	  // Calculate needed shift and apply it.
	  double c1  = 0.5 * ( (1. - r1 + r2) * l34 / l12 - (1. - r3 + r4) );
	  double c2  = 0.5 * ( (1. + r1 - r2) * l34 / l12 - (1. + r3 - r4) );
	  Vec4   pSh = c1 * p1Move - c2 * p2Move;
	  p1Move    += pSh;
	  p2Move    -= pSh;
	  return true;
	}

	//--------------------------------------------------------------------------

	// Create two vectors that are perpendicular to both input vectors.

	pair<Vec4,Vec4> getTwoPerpendicular(const Vec4& v1, const Vec4& v2) {

	  // One perpendicular vector from three-dimensional cross-product.
	  Vec4 nPerp( cross3(v1,v2) );
	  double TINY = std::numeric_limits<double>::epsilon();
	  if ( abs(nPerp.pAbs()) < TINY) {
		Vec4 aux;
		if (v1.x() != 0.)      aux.p(v1.yy,v1.x(),v1.z(),v1.t());
		else if (v1.y() != 0.) aux.p(v1.x(),v1.z(),v1.y(),v1.t());
		else if (v1.z() != 0.) aux.p(v1.z(),v1.y(),v1.x(),v1.t());
		nPerp.p( cross3(v1,aux) );
	  }
	  nPerp /= abs(nPerp.pAbs());

	  // Second perpendicular vector from four-dimensional cross-product.
	  Vec4 lPerp( cross4(v1,v2,nPerp) );
	  lPerp /= sqrt(abs(lPerp.m2Calc()));
	  return make_pair(nPerp,lPerp);
	}

}

//End of file
