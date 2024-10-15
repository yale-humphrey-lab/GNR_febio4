    #include "Elastin.h"
#include "Conf.h"
#include <fstream>

Elastin::Elastin():Constituent(), m_spatial_modulus(0.0)
{
}

Elastin::Elastin(double mass, double g_circ, double g_axial, double modulus, double Jdep, double bulkLM)
	: Constituent(mass, modulus, g_circ, false, false), m_Jdep(Jdep), m_bulkLM(bulkLM), m_G_axial(g_axial), m_spatial_modulus(0.0)
{
}
Elastin::Elastin(const Elastin& other) : Constituent(other.m_mass_frac, other.m_modulus, other.m_G, false, 0.0),
m_G_cric(other.m_G_cric), m_G_axial(other.m_G_axial), m_Jdep(other.m_Jdep), m_bulkLM(other.m_bulkLM)
{
}
void Elastin::calcEnergy(Conf* pCnf)
{
	mat3d ip;
	switch (pCnf->getMode())
	{
	case MODE_PRE:
		ip = pCnf->getF() * m_G_;
	case MODE_PROCESS:
		ip = pCnf->getF() * m_G_; //Todo
	case MODE_POST:
		ip = pCnf->getF() * m_G_;
	}
	double energy = m_mass_frac * (m_modulus / 2.0 * ((ip).dotdot(ip) - 3.0)); //decide if the mass frac is in this level or conf level
}
bool Elastin::isPenalty()
{
	return true;
}
mat3ds Elastin::calcActiveStress(Conf* pCnf)
{
	return mat3ds(0.0);
}
mat3ds Elastin::calcPassiveStress(Conf* pCnf)
{

	mat3ds Se(0.0);
	switch (pCnf->getMode())
	{
	case MODE_PRE:
		Se = (m_modulus * m_G_ * m_G_).sym();
	case MODE_PROCESS:
		Se = (m_modulus * m_G_ * m_G_).sym(); //Todo
	case MODE_POST:
		Se = (m_modulus * m_G_ * m_G_).sym();
	}
	m_psv_cauchy_stress = Se;
	return Se;
}
tens4dmm Elastin::calcTangentPassiveComponent(Conf* pCnf)
{
	//tens4ds ce(0.0);
	//return m_spatial_modulus;

	return tens4dmm(m_spatial_modulus);
}
tens4dmm Elastin::calcTangentActiveComponent(Conf* pCnf)
{
	return tens4dmm(0.0);
}
double Elastin::calcScalarPenalty(Conf* pCnf)
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
tens4dmm Elastin::calcTangentPenalty(Conf* pCnf)
{
	const double J = pCnf->getJ();
	const double lm = m_bulkLM * m_modulus;

	const mat3dd  I(1.0);
	const tens4ds IxI = dyad1s(I);
	const tens4ds IoI = dyad4s(I);
	
	tens4ds pn(0.0);

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
void Elastin::init_G(Conf* pCnf)
{
	vec3d N0 = pCnf->getCurrentLayer()->getBaseVector()[0];
	vec3d N1 = pCnf->getCurrentLayer()->getBaseVector()[1];
	vec3d N2 = pCnf->getCurrentLayer()->getBaseVector()[2];
	
	m_G_ = 1.0 / m_G / m_G_axial * dyad(N0)   + m_G * dyad(N1) +  m_G_axial * dyad(N2);

}
void Elastin::init(Conf* pCnf)
{
	this->init_G(pCnf);
}