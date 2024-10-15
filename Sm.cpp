#include "Sm.h"
#include "Conf.h"
#include <fstream>

Sm::Sm(): Constituent(), m_lr(0.0), m_lt(0.0), m_lto(0.0), m_lmt2(0.0), m_Cratio(0.0)
{
}
Sm::Sm(double mass, double G, double modulus_a, double modulus_b, double tmax, double CB, double max_c, double min_c)
	: Constituent(mass, modulus_a, G, true, false), m_modulus_b(modulus_b), m_Tmax(tmax), m_CB(CB), m_CS(0.5 * CB), m_max_contraction(max_c), m_min_contraction(min_c)
{
	m_lr, m_lt, m_lto, m_lmt2, m_Cratio = 0.0;
	
}
Sm::Sm(double mass, double G, double modulus_a, double modulus_b, double tmax, double CB, double CS, double max_c, double min_c)
	: Constituent(mass, modulus_a, G, true, false), m_modulus_b(modulus_b), m_Tmax(tmax), m_CB(CB), m_CS(CS), m_max_contraction(max_c), m_min_contraction(min_c)
{
	m_lr, m_lt, m_lto, m_lmt2, m_Cratio = 0.0;
}
Sm::Sm(const Sm& other) : Constituent(other.m_mass_frac, other.m_modulus, other.m_G, true, 0.0),
m_modulus_b(other.m_modulus_b), m_Tmax(other.m_Tmax), m_CB(other.m_CB), m_CS(other.m_CS), m_max_contraction(other.m_max_contraction), m_min_contraction(other.m_min_contraction)
{
}
void Sm::init(Conf* pCnf)
{
	Mode  m = pCnf->getMode();
	m_fib_N = pCnf->getCurrentLayer()->getBaseVector()[1];

	
	if (m == MODE_PRE)// || m == MODE_PROCESS)
	{	
		m_lto = (pCnf->getF() * m_fib_N).norm();
		m_lt = m_lto;
		m_tent = dyad(pCnf->getF() * m_fib_N);
	}
	else if (m == MODE_PROCESS)
	{
		m_lt = (pCnf->getF() * (pCnf->getRefConf()->getFinv() * m_fib_N)).norm();
		m_lto = (pCnf->getRefConf()->getF() * m_fib_N).norm();

		m_lr = (pCnf->getF() * (pCnf->getRefConf()->getFinv() * pCnf->getCurrentLayer()->getBaseVector()[0])).norm();
		m_tent = dyad(pCnf->getF() * (pCnf->getRefConf()->getFinv() * m_fib_N));
	}
	else //if (m == MODE_POST)
	{
		mat3ds Uo; mat3d Ro; (pCnf->getRefConf()->getF()).right_polar(Ro, Uo);	// Uo from polar decomposition
		mat3ds Uh; mat3d Rh; (pCnf->getHomConf()->getF()).right_polar(Rh, Uh);	// Uh from polar decomposition

		m_lto = (pCnf->getF() * (Uh.inverse() * (Uo * m_fib_N))).norm();
		m_tent = dyad(pCnf->getF() * (Uh.inverse() * (Uo * m_fib_N)));
		m_lr = (pCnf->getHomConf()->getF() * (Uh.inverse() * (Uo * pCnf->getCurrentLayer()->getBaseVector()[0]))).norm();
	}	
	m_lmt2 = (m_G * m_lto) * (m_G * m_lto);
}
void Sm::calcEnergy(Conf* pCnf)
{
	double energy = m_modulus * (m_modulus / (4.0 * m_modulus_b) * (exp(m_modulus_b * (m_lmt2 - 1.0) * (m_lmt2 - 1.0)) - 1.0));
}
bool Sm::isActiveStress()
{
	return true;
}
mat3ds Sm::calcActiveStress(Conf* pCnf)
{
	mat3ds Sa(0.0);
	Sa = m_Tmax * (1.0 - exp(-m_CB * m_CB)) * (1.0 - pow((m_max_contraction - 1.0) / (m_max_contraction - m_min_contraction), 2)) * (m_lto * m_lto) * dyad(m_fib_N);
	
	if (pCnf->getMode() == MODE_PROCESS)
	{
		mat3d  uo(pCnf->getRefConf()->getU());
		m_active_cauchy_no_decay = 1.0 / pCnf->getRefConf()->getJ() * (uo * (Sa * uo)).sym();
		m_act_cauchy_stress = m_active_cauchy_no_decay;
		//applying decay
		const double rr = pCnf->getCurrentLayer()->getRRatio();
		double rIrIo = rr * m_lt - (rr - 1) * m_lr;				
		double Cratio = m_CB - m_CS * (pow(rIrIo, -3) - 1.0);

		if (Cratio > 0)
		{
			m_act_cauchy_stress = (1.0 - exp(-Cratio * Cratio)) / (1.0 - exp(-m_CB * m_CB)) * m_active_cauchy_no_decay;
		}
		else
		{
			m_act_cauchy_stress.zero();
		}
		
	}
	return Sa;
}
double Sm::calcScalarPenalty(Conf* pCnf)
{
	return 0.0;
}
tens4dmm Sm::calcTangentPenalty(Conf* pCnf)
{
	return tens4dmm(0.0);
}
mat3ds Sm::calcPassiveStress(Conf* pCnf)
{
	vec3d v;
	mat3ds u;
	double j_;
	double Jhom_to_j = 1;
	switch (pCnf->getMode())
	{
	case MODE_PRE:
		v = m_fib_N;
		j_ = pCnf->getJ();
		u = pCnf->getU();
		break;
	case MODE_PROCESS:
		v = m_fib_N;
		j_ = pCnf->getRefConf()->getJ();
		u = pCnf->getRefConf()->getU();
		break;
	case MODE_POST:
		u = pCnf->getRefConf()->getU();
		v = pCnf->getHomConf()->getUinv() * (u * m_fib_N);
		j_ = pCnf->getRefConf()->getJ();
		Jhom_to_j = pCnf->getHomConf()->getJ() / j_;
		
		break;
	}
	const mat3ds Sm = Jhom_to_j * (m_modulus * (m_lmt2 - 1.0) * exp(m_modulus_b * (m_lmt2 - 1.0) * (m_lmt2 - 1.0)) * (m_G * m_G) * dyad(v));

	mat3d u_(u);
	m_psv_cauchy_stress = 1.0 / j_ * (u_ * (Sm * u_)).sym();

	return Sm;
}
/*Tangent*/
tens4dmm Sm::calcTangentPassiveComponent(Conf* pCnf)
{
	
	double Jhom_to_j = 1;

	const double exp_term = m_modulus_b * (m_lmt2 - 1.0) * (m_lmt2 - 1.0);

	if (pCnf->getMode() == MODE_POST)
	{
		Jhom_to_j = pCnf->getHomConf()->getJ() / pCnf->getRefConf()->getJ();
	}

	tens4ds tSm = Jhom_to_j * (2.0 * m_modulus * (1.0 + 2.0 * exp_term) * exp(exp_term) * pow(m_G, 4) * dyad1s(m_tent));

	return tens4dmm(tSm);
}
tens4dmm Sm::calcTangentActiveComponent(Conf* pCnf)
{
	Mode m = pCnf->getMode();
	if (m == MODE_PRE)
	{	
		tens4ds tCa = 2.0 * m_Tmax * (1.0 - exp(-m_CB * m_CB)) * (1.0 - pow((m_max_contraction - 1.0) / (m_max_contraction - m_min_contraction), 2)) * dyad1s(m_tent);
		return tens4dmm(tCa);
	}
	if (m == MODE_PROCESS)
	{
		mat3d F = pCnf->getF();
		mat3d R = pCnf->getR();
		mat3d Fio = pCnf->getRefConf()->getFinv();
		vec3d N[3]; pCnf->getCurrentLayer()->getBaseVector(N);
		mat3ds tenr = dyad(F * (Fio * N[0])); // const						// Fio needed for consistency (from computation of lr)
		mat3ds tent = dyad(F * (Fio * N[1])); // const
		mat3ds tenz = dyad(F * (Fio * N[2])); // const

		const double rr = pCnf->getCurrentLayer()->getRRatio();
		double rIrIo = rr * m_lt - (rr - 1) * m_lr;				// rIrIo -> rIorIo = 1 for F -> Fo
		double Cratio = m_CB - m_CS * (pow(rIrIo, -3) - 1.0); // const+
		// Contribution due to the ratio of vasocontrictors to vasodilators in the active stress
		tens4dmm saoxnrr = dyad1mm((R * m_active_cauchy_no_decay * R.transpose()).sym(), tenr); // const
		tens4dmm saoxntt = dyad1mm((R * m_active_cauchy_no_decay * R.transpose()).sym(), tent); // const

		// 1/J * FoF : [ J * phim * 1/(1.0-exp(-CB*CB)) * (Ui*sao*Ui) x d(1-exp(-Cratio^2))/d(C/2) ] : (Ft)o(Ft)
		tens4dmm cass = 6.0 * Cratio * m_CS * pow(rIrIo, -4) * exp(-Cratio * Cratio) /
			(1.0 - exp(-m_CB * m_CB)) * (rr / m_lt * saoxntt - (rr - 1) / m_lr * saoxnrr); // const

		return (cass);
	}
	if (m == MODE_POST)
	{
		return tens4dmm(0.0);
	}	
}
