#include "Collagen.h"
#include "Conf.h"
#include <fstream>


Collagen::Collagen():Constituent(), m_lt(0.0), m_lz(0.0), m_lp(0.0), m_ln(0.0), m_lct2(0.0), m_lcz2(0.0), m_lcp2(0.0), m_lcn2(0.0),
m_fib_N1(vec3d(0.0)), m_fib_N2(vec3d(0.0)), m_fib_Np(vec3d(0.0)), m_fib_Nn(vec3d(0.0)){}

Collagen::Collagen(double mass, double G, double modulus_a, double modulus_b, double c_frc, double a_frc, double angle, double reor)
	: Constituent(mass, modulus_a, G, true, angle, reor), m_modulus_b(modulus_b), 
	m_lt(0.0), m_lz(0.0), m_lp(0.0), m_ln(0.0), m_lct2(0.0), m_lcz2(0.0), m_lcp2(0.0), m_lcn2(0.0),
	m_fib_N1(vec3d(0.0)), m_fib_N2(vec3d(0.0)), m_fib_Np(vec3d(0.0)), m_fib_Nn(vec3d(0.0)),
	m_circ_frac(c_frc), m_axial_frac(a_frc), m_diag_frac(0.5 * (1.0 - c_frc - a_frc)){}

Collagen::Collagen(const Collagen& other) : Constituent(other.m_mass_frac, other.m_modulus, other.m_G, true, other.m_angle, other.m_reorientation),
 m_modulus_b(other.m_modulus_b),m_circ_frac(other.m_circ_frac), m_axial_frac(other.m_axial_frac), m_diag_frac(other.m_diag_frac),
	m_lt(other.m_lt), m_lz(other.m_lz), m_lp(other.m_lp), m_ln(other.m_ln), m_lct2(other.m_lct2), m_lcz2(other.m_lcz2), m_lcp2(other.m_lcp2), m_lcn2(other.m_lcn2),
	m_fib_N1(other.m_fib_N1), m_fib_N2(other.m_fib_N2), m_fib_Np(other.m_fib_Np), m_fib_Nn(other.m_fib_Nn)
{
}
double Collagen::calcScalarPenalty(Conf* pCnf)
{
	return 0.0;
}
tens4dmm Collagen::calcTangentPenalty(Conf* pCnf)
{
	return tens4dmm(0.0);
}
mat3ds Collagen::calcActiveStress(Conf* pCnf)
{
	return mat3ds(0.0);
}
mat3ds Collagen::calcPassiveStress(Conf* pCnf)
{
//	init(pCnf);

	double Jhmo_Jref_ratio = 1;
	
	Mode m = pCnf->getMode();

	if (m == MODE_POST)
	{
		mat3ds Uo = pCnf->getRefConf()->getU();
		mat3ds Uh_inv = pCnf->getHomConf()->getUinv();
		Jhmo_Jref_ratio = pCnf->getHomConf()->getJ() / pCnf->getRefConf()->getJ();

		const mat3ds Sc = Jhmo_Jref_ratio *
			(	m_modulus * (m_lct2 - 1.0) * exp(m_modulus_b * (m_lct2 - 1.0) * (m_lct2 - 1.0)) * (m_G * m_G) * dyad(Uh_inv * (Uo * m_fib_N1)) * m_circ_frac +
				m_modulus * (m_lcz2 - 1.0) * exp(m_modulus_b * (m_lcz2 - 1.0) * (m_lcz2 - 1.0)) * (m_G * m_G) * dyad(Uh_inv * (Uo * m_fib_N2)) * m_axial_frac +
				m_modulus * (m_lcp2 - 1.0) * exp(m_modulus_b * (m_lcp2 - 1.0) * (m_lcp2 - 1.0)) * (m_G * m_G) * dyad(Uh_inv * (Uo * m_fib_Np)) * m_diag_frac +
				m_modulus * (m_lcn2 - 1.0) * exp(m_modulus_b * (m_lcn2 - 1.0) * (m_lcn2 - 1.0)) * (m_G * m_G) * dyad(Uh_inv * (Uo * m_fib_Nn)) * m_diag_frac);

		return Sc;
	}
	
	const mat3ds Sc = Jhmo_Jref_ratio * 
		(m_modulus * (m_lct2 - 1.0) * exp(m_modulus_b * (m_lct2 - 1.0) * (m_lct2 - 1.0)) * (m_G * m_G) * dyad(m_fib_N1) * m_circ_frac +
		 m_modulus * (m_lcz2 - 1.0) * exp(m_modulus_b * (m_lcz2 - 1.0) * (m_lcz2 - 1.0)) * (m_G * m_G) * dyad(m_fib_N2) * m_axial_frac +
		 m_modulus * (m_lcp2 - 1.0) * exp(m_modulus_b * (m_lcp2 - 1.0) * (m_lcp2 - 1.0)) * (m_G * m_G) * dyad(m_fib_Np) * m_diag_frac +
		 m_modulus * (m_lcn2 - 1.0) * exp(m_modulus_b * (m_lcn2 - 1.0) * (m_lcn2 - 1.0)) * (m_G * m_G) * dyad(m_fib_Nn) * m_diag_frac);

	// Calc Cauchy stress
	if (m == MODE_PRE)
	{
		mat3d u(pCnf->getU());
		m_psv_cauchy_stress = 1.0 / pCnf->getJ() * (u * (Sc * u)).sym();
	}
	else if (m == MODE_PROCESS)
	{
		mat3d u(pCnf->getRefConf()->getU());
		m_psv_cauchy_stress = 1.0 / pCnf->getRefConf()->getJ() * (u * (Sc * u)).sym();
	}
	else // MODE_POST: currently not in use
	{
//		mat3d u(pCnf->getPrvConf()->getU());
//		m_psv_cauchy_stress = 1.0 / pCnf->getPrvConf()->getJ() * (u * (Sc * u)).sym();
		
	}
	return Sc;
}
tens4dmm Collagen::calcTangentPassiveComponent(Conf* pCnf)
{
	const double exp_term1 = m_modulus_b * (m_lct2 - 1.0) * (m_lct2 - 1.0);
	const double exp_term2 = m_modulus_b * (m_lcz2 - 1.0) * (m_lcz2 - 1.0);
	const double exp_term3 = m_modulus_b * (m_lcp2 - 1.0) * (m_lcp2 - 1.0);
	const double exp_term4 = m_modulus_b * (m_lcn2 - 1.0) * (m_lcn2 - 1.0);

	double jHom_to_J = 1;

	switch (pCnf->getMode())
	{
	case MODE_PRE:
		break;
	case MODE_PROCESS:
		break;
	case MODE_POST: 
		jHom_to_J = pCnf->getHomConf()->getJ() / pCnf->getRefConf()->getJ();
		break;
	} 
	tens4ds tCf = jHom_to_J * 
		(2.0 * m_modulus * (1.0 + 2.0 * exp_term1) * exp(exp_term1) * pow(m_G, 4) * dyad1s(m_tent) * m_circ_frac +
		2.0 * m_modulus * (1.0 + 2.0 * exp_term2) * exp(exp_term2) * pow(m_G, 4) * dyad1s(m_tenz) * m_axial_frac +
		2.0 * m_modulus * (1.0 + 2.0 * exp_term3) * exp(exp_term3) * pow(m_G, 4) * dyad1s(m_tenp) * m_diag_frac +
		2.0 * m_modulus * (1.0 + 2.0 * exp_term4) * exp(exp_term4) * pow(m_G, 4) * dyad1s(m_tenn) * m_diag_frac);

	return tens4dmm(tCf);
}
tens4dmm Collagen::calcTangentActiveComponent(Conf* pCnf)
{
	return tens4dmm(0.0);
}
tens4dmm Collagen::calcTangetOrientationComponent(Conf* pCnf)
{
	const double scphato = m_modulus * (m_lcp2 - 1.0) * exp(m_modulus_b * (m_lcp2 - 1.0) * (m_lcp2 - 1.0)) * (m_G * m_G);	// Constant stress magnitude at constituent level
	const double scnhato = m_modulus * (m_lcn2 - 1.0) * exp(m_modulus_b * (m_lcn2 - 1.0) * (m_lcn2 - 1.0)) * (m_G * m_G);
	
	const vec3d dNpdta = (m_fib_N1 - m_fib_N2 * tan(m_angle)) * pow(1 + pow(tan(m_angle), 2), -1.5);	// d(Np)/d(tan(alpha))
	const vec3d dNndta = (m_fib_N1 + m_fib_N2 * tan(m_angle)) * pow(1 + pow(tan(m_angle), 2), -1.5);

	const mat3d R = pCnf->getR();
	const mat3ds Uo = pCnf->getRefConf()->getU();
	const double Jo = pCnf->getRefConf()->getJ();
	const double cur_lt = (pCnf->getF() * (pCnf->getRefConf()->getFinv() * m_fib_N1)).norm();
	const double cur_lz = (pCnf->getF() * (pCnf->getRefConf()->getFinv() * m_fib_N2)).norm();
	const mat3ds ten1 = 1.0 / Jo * dyads(R * (Uo * dNpdta), R * (Uo * m_fib_Np));
	const mat3ds ten2 = 1.0 / Jo * dyads(R * (Uo * dNndta), R * (Uo * m_fib_Nn));
	const mat3ds ten3 = m_reorientation * tan(m_angle) * (1.0 / (cur_lt * cur_lt) * m_tent - 1.0 / (cur_lz * cur_lz) * m_tenz);		// 2*d(tan(alpha))/d(C) : (Ft)o(Ft)
	
	tens4dmm cpnss(0.0);
	cpnss += (m_mass_frac * m_diag_frac) * scphato * dyad1mm(ten1, ten3);
	cpnss += (m_mass_frac * m_diag_frac) * scnhato * dyad1mm(ten2, ten3);
	
	return cpnss;
}
void Collagen::init(Conf* pCnf)
{
	m_fib_N1 = pCnf->getCurrentLayer()->getBaseVector()[1];
	m_fib_N2 = pCnf->getCurrentLayer()->getBaseVector()[2];
	m_fib_Np = m_fib_N1 * sin(m_angle) + m_fib_N2 * cos(m_angle);		// Original diagonal fiber direction
	m_fib_Nn = m_fib_N1 * sin(m_angle) - m_fib_N2 * cos(m_angle);		// idem for symmetric

	const mat3d F = pCnf->getF();

	Mode m = pCnf->getMode();

	if (m == MODE_PRE)
	{

		m_lt = (F * m_fib_N1).norm();		m_lz = (F * m_fib_N2).norm();
		m_lp = (F * m_fib_Np).norm();		m_ln = (F * m_fib_Nn).norm();

		m_tent = dyad(F * m_fib_N1);		m_tenz = dyad(F * m_fib_N2);
		m_tenp = dyad(F * m_fib_Np);		m_tenn = dyad(F * m_fib_Nn);
	}
	else if (m == MODE_POST)
	{
		mat3d Ro, Rh;
		mat3ds Uo, Uh;

		pCnf->getRefConf()->getF().right_polar(Ro, Uo);
		pCnf->getHomConf()->getF().right_polar(Rh, Uh);

		m_lt = (F * (Uh.inverse() * (Uo * m_fib_N1))).norm();		m_lz = (F * (Uh.inverse() * (Uo * m_fib_N2))).norm();

		m_lp = (F * (Uh.inverse() * (Uo * m_fib_Np))).norm();		m_ln = (F * (Uh.inverse() * (Uo * m_fib_Nn))).norm();

		m_tent = dyad(F * (Uh.inverse() * (Uo * m_fib_N1)));		m_tenz = dyad(F * (Uh.inverse() * (Uo * m_fib_N2)));
		m_tenp = dyad(F * (Uh.inverse() * (Uo * m_fib_Np)));		m_tenn = dyad(F * (Uh.inverse() * (Uo * m_fib_Nn)));
	}
	else // MODE_PROCESS
	{
		
		const mat3d F_ = pCnf->getRefConf()->getF();
		double org_angle = pCnf->getRefConf()->getCurrentLayer()->getTurnOverConstituents().at(0)->getAngle();
		//vec3d fib_N1 = pCnf->getRefConf()->getCurrentLayer()->getBaseVector()[1];
		//vec3d fib_N2 = pCnf->getRefConf()->getCurrentLayer()->getBaseVector()[2];
		vec3d fib_Np = m_fib_N1 * sin(org_angle) + m_fib_N2 * cos(org_angle);		// Original diagonal fiber direction
		vec3d fib_Nn = m_fib_N1 * sin(org_angle) - m_fib_N2 * cos(org_angle);		// idem for symmetric

		m_lt = (F_ * m_fib_N1).norm();
		m_lz = (F_ * m_fib_N2).norm();
		m_lp = (F_ * fib_Np).norm();
		m_ln = (F_ * fib_Nn).norm();

		m_tent = dyad(F * (pCnf->getRefConf()->getFinv() * m_fib_N1));		m_tenz = dyad(F * (pCnf->getRefConf()->getFinv() * m_fib_N2));
	}

	m_lct2 = (m_G * m_lt) * (m_G * m_lt);		m_lcz2 = (m_G * m_lz) * (m_G * m_lz);
	m_lcp2 = (m_G * m_lp) * (m_G * m_lp);		m_lcn2 = (m_G * m_ln) * (m_G * m_ln);
}
void Collagen::calcEnergy(Conf* pCnf)
{

	double energy = (m_modulus / (4.0 * m_modulus_b) * (exp(m_modulus_b * (m_lct2 - 1.0) * (m_lct2 - 1.0)) - 1.0) * m_circ_frac +
		m_modulus / (4.0 * m_modulus_b) * (exp(m_modulus_b * (m_lcz2 - 1.0) * (m_lcz2 - 1.0)) - 1.0) * m_axial_frac +
		m_modulus / (4.0 * m_modulus_b) * (exp(m_modulus_b * (m_lcp2 - 1.0) * (m_lcp2 - 1.0)) - 1.0) * m_diag_frac +
		m_modulus / (4.0 * m_modulus_b) * (exp(m_modulus_b * (m_lcn2 - 1.0) * (m_lcn2 - 1.0)) - 1.0) * m_diag_frac);
}
void Collagen::setMassFrac(double massFraction, double circMassFraction, double axialMassFraction)
{
	setMassFraction(massFraction);

	m_circ_frac = circMassFraction;
	m_axial_frac = axialMassFraction;

	m_diag_frac = 0.5 * (1.0 - m_circ_frac - m_axial_frac);
}
