#ifndef COLLAGEN_H
#define COLLAGEN_H

#include "Constituent.h"

class Collagen : public Constituent
{
public:
	Collagen();
	Collagen(double mass, double G, double modulus_a, double modulus_b, double c_frc, double a_frc, double angle, double reor = 0.0);
	Collagen(const Collagen& other);
	Collagen* clone() const override { return new Collagen(*this); }
	//////////////////////////////////////
	void init(Conf* pCnf) override;

	bool isTurnOver() override { return true; }
	bool isActiveStress() override { return false; }
	bool isPenalty() override { return false; }
	
	void calcEnergy(Conf* pCnf) override;

	double calcScalarPenalty(Conf* pCnf) override;
	tens4dmm calcTangentPenalty(Conf* pCnf) override;

	mat3ds calcActiveStress(Conf* pCnf) override;
	mat3ds calcPassiveStress(Conf* pCnf) override;
	tens4dmm calcTangentPassiveComponent(Conf* pCnf) override;
	tens4dmm calcTangentActiveComponent(Conf* pCnf) override;
	tens4dmm calcTangetOrientationComponent(Conf* pCnf) override;

	void setMassFrac(double massFraction, double circMassFraction, double axialMassFraction);

private:
	double m_modulus_b;

	vec3d  m_fib_N1, m_fib_N2, m_fib_Np, m_fib_Nn;
	double m_circ_frac, m_axial_frac, m_diag_frac;
	double m_lt, m_lz, m_lp, m_ln;
	double m_lct2, m_lcz2, m_lcp2, m_lcn2;
	mat3ds m_tent, m_tenz, m_tenp, m_tenn;

};

#endif