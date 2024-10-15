#include "Constituent.h"
#include "FECore/FEAnalysis.h"
#include "FECore/FEModel.h"	
#include "Conf.h"
#include <fstream>

Constituent::Constituent(): m_psv_cauchy_stress(mat3ds(0.0)), m_act_cauchy_stress(mat3ds(0.0)), m_lm(0.0)
{
}

Constituent::Constituent(double mass, double modulus, double G, bool to, double angle, double reorient) :
	m_modulus(modulus), m_mass_frac(mass), m_G(G), m_turn_over(to), m_angle(angle), m_reorientation(reorient),
	m_psv_cauchy_stress(mat3ds(0.0)), m_act_cauchy_stress(mat3ds(0.0)), m_lm(0.0)
{
}
Constituent::Constituent(double mass, double modulus, double G, bool to, double reorient) :
	m_modulus(modulus), m_mass_frac(mass), m_G(G), m_turn_over(to), m_reorientation(reorient), 
	m_psv_cauchy_stress(mat3ds(0.0)) ,m_act_cauchy_stress(mat3ds(0.0)), m_lm(0.0), m_angle(0.0)
{
	// Initialize other members here
}
Constituent::~Constituent()
{
}
void Constituent::setParameters(double mass, double modulus, double G, bool to, double reorient)
{
	m_mass_frac = mass;
	m_modulus = modulus;
	m_G = G;
	m_turn_over = to;
	m_reorientation = reorient;
}
void Constituent::setParameters(double mass, double modulus, double G, bool to, double angle, double reorient)
{
	m_mass_frac = mass;
	m_modulus = modulus;
	m_G = G;
	m_turn_over = to;
	m_angle = angle;
	m_reorientation = reorient;
}
//const double Constituent::getEnergy()
//{
//	return m_energy;
//}
void Constituent::recalcAngle(Conf* pCnf)
{
	if (pCnf->getMode() == MODE_PROCESS)
	{
		//const double lt = (pCnf->getF() * (pCnf->getPrvConf()->getFinv() * pCnf->getCurrentLayer()->getBaseVector()[1])).norm();
		//const double lz = (pCnf->getF() * (pCnf->getPrvConf()->getFinv() * pCnf->getCurrentLayer()->getBaseVector()[2])).norm();
		const double lt = (pCnf->getF() * (pCnf->getRefConf()->getFinv() * pCnf->getCurrentLayer()->getBaseVector()[1])).norm();
		const double lz = (pCnf->getF() * (pCnf->getRefConf()->getFinv() * pCnf->getCurrentLayer()->getBaseVector()[2])).norm();
		const double ore = pCnf->getRefConf()->getCurrentLayer()->getTurnOverConstituents().at(0)->getReor();
		double n_angle = atan(tan(m_angle) * pow(lt / lz, ore));
		m_angle = n_angle;
	}
		if (pCnf->getMode() == MODE_POST)
		{
			const double lt = (pCnf->getHomConf()->getF() * (pCnf->getPrvConf()->getFinv() * pCnf->getCurrentLayer()->getBaseVector()[1])).norm();
			const double lz = (pCnf->getHomConf()->getF() * (pCnf->getPrvConf()->getFinv() * pCnf->getCurrentLayer()->getBaseVector()[2])).norm();
			const double ore = pCnf->getRefConf()->getCurrentLayer()->getTurnOverConstituents().at(0)->getReor();
			m_angle = atan(tan(m_angle) * pow(lt / lz, ore));
		}
}
double Constituent::getFraction()
{
	return m_mass_frac;
}
void Constituent::setMassFraction(double frac)
{
	m_mass_frac = frac;
}
void Constituent::setAngle(double alpha)
{
	m_angle = alpha;
}