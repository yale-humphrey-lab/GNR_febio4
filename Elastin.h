#ifndef ELASTIN_H
#define ELASTIN_H

#include "Constituent.h"
#define HOM_JDEP 0.9999
#include <fstream>

class Elastin : public Constituent
{
public:
	Elastin();
	Elastin(double mass, double g_circ, double g_axial, double modulus, double Jdep, double bulkLM);
	Elastin(const Elastin& other);
	Elastin* clone() const override { return new Elastin(*this); }
	void init(Conf* pCnf) override;
	
	void calcEnergy(Conf* pCnf) override;
	
	bool isPenalty() override;
	bool isTurnOver() override { return false; }
	bool isActiveStress() override { return false; }

	mat3ds calcActiveStress(Conf* pCnf) override; 	
	mat3ds calcPassiveStress(Conf* pCnf) override;
	tens4dmm calcTangentPassiveComponent(Conf* pCnf) override;
	tens4dmm calcTangentActiveComponent(Conf* pCnf) override;
	tens4dmm calcTangetOrientationComponent(Conf* pCnf) override { return tens4dmm(0.0); }

	double calcScalarPenalty(Conf* pCnf) override;
	tens4dmm calcTangentPenalty(Conf* pCnf) override;

	void init_G(Conf* pCnf);
	void setJdep(double Jp) { m_Jdep = HOM_JDEP; }
	const tens4ds getSpatialModulus() { return m_spatial_modulus; }
	
	void print();

private:
	double m_Jdep;        // Near incompressibility
	double m_bulkLM;      // Incompressibility multiplier

	const tens4ds m_spatial_modulus; //ce
	mat3ds m_G_;
	double m_G_cric;
	double m_G_axial;
};

#endif