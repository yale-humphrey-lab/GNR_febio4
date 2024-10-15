#include "GaGs.h"
#include "Conf.h"
#include <fstream>
Gags::Gags() :Constituent(), m_spatial_modulus(0.0)
{
}

Gags::Gags(double mass, double g_radial , double modulus, double Jdep, double bulkLM)
	: Constituent(mass, modulus, g_radial, true, false), m_Jdep(Jdep), m_bulkLM(bulkLM), m_spatial_modulus(0.0)
{
	m_psv_cauchy_stress.zero();
}
Gags::Gags(const Gags& other) : Constituent(other.m_mass_frac, other.m_modulus, other.m_G, true, 0.0),
m_G_radial(other.m_G_radial), m_Jdep(other.m_Jdep), m_bulkLM(other.m_bulkLM)
{
}
void Gags::calcEnergy(Conf* pCnf)
{
	mat3d ip;
	switch (pCnf->getMode())
	{
	case MODE_PRE:
		ip = pCnf->getF() * m_G_;
	case MODE_PROCESS:
		ip = pCnf->getF() * m_G_;
	case MODE_POST:
		ip = pCnf->getF() * m_G_;
	}
	double energy = m_mass_frac * (m_modulus / 2.0 * ((ip).dotdot(ip) - 3.0)); //decide if the mass frac is in this level or conf level
}
bool Gags::isPenalty()
{
	return true;
}
mat3ds Gags::calcActiveStress(Conf* pCnf)
{
	return mat3ds(0.0);
}
mat3ds Gags::calcPassiveStress(Conf* pCnf)
{
	mat3ds Se(0.0);

	Se = (m_modulus * m_G_ * m_G_).sym();
	m_psv_cauchy_stress = Se;
	return Se;
}
tens4dmm Gags::calcTangentPassiveComponent(Conf* pCnf)
{
	return tens4dmm(m_spatial_modulus);
}
tens4dmm Gags::calcTangentActiveComponent(Conf* pCnf)
{
	return tens4dmm(0.0);
}
double Gags::calcScalarPenalty(Conf* pCnf)
{
	switch (pCnf->getMode())
	{
	case MODE_PRE:
		return m_bulkLM * m_modulus * log(m_Jdep * pCnf->getJ());
	case MODE_PROCESS:
		return 0.0;
	case MODE_POST:
		return m_bulkLM * m_modulus * log(0.9999 * pCnf->getJ() / pCnf->getHomConf()->getJ());
	}
}
tens4dmm Gags::calcTangentPenalty(Conf* pCnf)
{
	const double J = pCnf->getJ();
	const double lm = m_bulkLM * m_modulus;

	const mat3dd  I(1.0);
	const tens4ds IxI = dyad1s(I);
	const tens4ds IoI = dyad4s(I);

	tens4ds pn(0.0);
	return tens4dmm(pn);

	switch (pCnf->getMode())
	{
	case MODE_PRE:
		pn = lm / J * (IxI - 2.0 * log(m_Jdep * J) * IoI);
		break;
	case MODE_PROCESS:
		//Todo
		break;
	case MODE_POST:
		double Jdep = 0.9999;
		pn = lm / J * (IxI - 2.0 * log(Jdep * J / pCnf->getHomConf()->getJ()) * IoI);
	}
	return tens4dmm(pn);
}
void Gags::init_G(Conf* pCnf)
{
	vec3d N0 = pCnf->getCurrentLayer()->getBaseVector()[0];
	vec3d N1 = pCnf->getCurrentLayer()->getBaseVector()[1];
	vec3d N2 = pCnf->getCurrentLayer()->getBaseVector()[2];

	//m_G_ = 1.0 / m_G / m_G_axial * dyad(N0) + m_G * dyad(N1) + m_G_axial * dyad(N2);
	m_G_ = m_G * dyad(N0) + 1/sqrt(m_G) * dyad(N1) + 1/sqrt(m_G) * dyad(N2);
}
void Gags::init(Conf* pCnf)
{
	this->init_G(pCnf);
}