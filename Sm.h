#ifndef _SM_H
#define _SM_H

#include "Constituent.h"

class Sm: public Constituent
{
public:
	Sm();
	Sm(double mass, double G, double modulus_a, double modulus_b, double tmax, double CB, double max_c, double min_c);
	Sm(double mass, double G, double modulus_a, double modulus_b, double tmax, double CB, double CS, double max_c, double min_c);
	Sm(const Sm& other);
	Sm* clone() const override { return new Sm(*this); }

	void init(Conf* pCnf) override;
	void calcEnergy(Conf* pCnf) override;

	bool isActiveStress() override;
	bool isTurnOver() override { return true; }
	bool isPenalty() override { return false; }

	void setActiveParam(double m) { m_Tmax = m; }
	double getActiveParam() { return m_Tmax; }

	double calcScalarPenalty(Conf* pCnf) override;
	tens4dmm calcTangentPenalty(Conf* pCnf) override;

	mat3ds calcPassiveStress(Conf* pCnf) override;
	mat3ds calcActiveStress(Conf* pCnf) override;

	tens4dmm calcTangentPassiveComponent(Conf* pCnf) override;
	tens4dmm calcTangentActiveComponent(Conf* pCnf) override;
	tens4dmm calcTangetOrientationComponent(Conf* pCnf) override { return tens4dmm(0.0); }

private:
	// Active
	double m_Tmax;      // Muscle tone, baseline
	double m_CB;         // Shape parameter s.t. (1-exp(-CB^2)) = 0.5
	double m_CS;         // Shape parameter s.t. (1-exp(-C^2)) = 0.0 for lt = 1/(1+CB/CS)^(1/3) = 0.7 and (1-exp(-C^2)) = 0.75 for lt = 2.0
	double m_max_contraction;       // Maximal contraction stretch
	double m_min_contraction;       // Minimal contraction stretch

	//& m_active_stress;

	double m_modulus_b;

	vec3d m_fib_N;
	mat3ds m_tent;

	double m_lr, m_lt, m_lto, m_lmt2;
	double m_Cratio;
	
	mat3ds m_active_cauchy_no_decay;

};

#endif