#include "Conf.h"
#include "Constituent.h"
#include "Insult.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <fstream>
#include <omp.h>
///////////////////
/*
		1. The stress of the non turn-over constituents is multiply by their original mass fraction and not the current.
		2. The actual mass fraction of each constituent is determined by collagen. Currently, it is hardcoded that the collegen is the first turnover constituent.
*/

Conf::Conf() : m_F(mat3d(0.0)), m_C(mat3ds(0.0)), m_J(1.0), m_current_t(0.0),
m_Ci(mat3ds(0.0)), m_Fi(mat3d(0.0)), m_U(mat3ds(0.0)), m_Ui(mat3ds(0.0)), m_R(mat3d(0.0)),
m_volumetric_stress(0.0),
m_conf_ref(nullptr), m_conf_hom(nullptr), m_conf_prv(nullptr),
m_currentLayer(nullptr), m_partial_t(11.0), m_mode(MODE_PRE)
{
	m_layers.clear();
}
Conf::Conf(const Conf& other): m_F(other.m_F), m_C(other.m_C), m_J(other.m_J), m_current_t(other.m_current_t),
		m_Ci(other.m_Ci), m_Fi(other.m_F.inverse()), m_U(other.m_U), m_Ui(other.m_Ui), m_R(other.m_R),
		m_volumetric_stress(other.m_volumetric_stress),
		m_conf_ref(nullptr), m_conf_hom(nullptr), m_conf_prv(nullptr),
		m_currentLayer(nullptr), m_partial_t(other.m_partial_t), m_mode(other.m_mode)
{
	m_layers.clear();
	if (other.m_layers.empty())
	{

		return;
	}
	for (std::size_t layerIndex = 0; layerIndex < other.m_layers.size(); layerIndex++)
	{
		m_layers[layerIndex] = new Layer(*other.m_layers[layerIndex]);
	}
}
Conf::~Conf()
{
	for (Layer* layer : m_layers)
	{
		delete layer;
	}
	m_layers.clear();
}
void Conf::setMode()
{
	if (m_current_t <= (1 + std::numeric_limits<double>::epsilon()))
	{
		m_mode = MODE_PRE;
	}
	else if (m_current_t <= (m_partial_t + std::numeric_limits<double>::epsilon()))
	{
		m_mode = MODE_PROCESS;
	}
	else
	{
		m_mode = MODE_POST;
	}
}
const Mode Conf::getMode()
{
	return m_mode;
}
void Conf::init(mat3d& F, mat3ds& C, double J, double t)
{
	m_F = F;					m_Fi = F.inverse();
	m_C = C;					m_Ci = C.inverse();

	m_F.right_polar(m_R, m_U);	m_Ui = m_U.inverse();

	m_J = J;

	//Times members
	m_current_t = t;
	setMode();
}
int Conf::setIrIo()
{
	try
	{
		const double clr = (m_F * (this->getRefConf()->getFinv() * getCurrentLayer()->getBaseVector()[0])).norm();	
		const double clt = (m_F * (this->getRefConf()->getFinv() * getCurrentLayer()->getBaseVector()[1])).norm();	
		getCurrentLayer()->setIrIo(clt, clr);
		return 1;
	}
	catch (...)
	{
		return 0;
	}
}
int Conf::setPosition(vec3d& v, bool tor)
{
	try
	{
		Layer* cL = getCurrentLayer();
		cL->setPosition(v);
		

		return 1;
	}
	catch (...)
	{

		return 0;
	}
}
void Conf::addLayer(Layer* input_layer)
{
	
	if (input_layer)
	{	
		m_layers.push_back(input_layer);
	}
}
void Conf::setCurrentLayer(Layer* input_layer)
{
	if (input_layer == nullptr)
	{
		std::string Strrsm = "D:/noG/InputLayerNullPtr";
		std::ofstream rf2(Strrsm);
		return;
	}
	m_currentLayer = input_layer;
}
Layer* Conf::getCurrentLayer()
{
	return m_currentLayer;
}
Layer* Conf::getCurrentLayer(vec3d& v)
{
	return m_layers[0];
}
Layer* Conf::getLayerByIndex(std::size_t layerIndex)
{
	return (m_layers)[layerIndex];
}

void Conf::setRefConf(Conf* pOther)
{
	m_conf_ref = pOther;
}
void Conf::setHomConf(Conf* pOther)
{

	m_conf_hom = pOther;
}
void Conf::setPreviousConf(Conf* pOther)
{
	m_conf_prv = pOther;

}

double Conf::recalcFibersAngle(Layer* curLayer)
{
	mat3d F(0.0);
	double orig_angle = 0.0;

	switch (m_mode)
	{
	case MODE_POST:
		F = m_conf_hom->getF();
		orig_angle = getHomConf()->getCurrentLayer()->getTurnOverConstituents().at(0)->getAngle();
		break;
	case MODE_PROCESS:
		F = getF();
		orig_angle = getRefConf()->getCurrentLayer()->getTurnOverConstituents().at(0)->getAngle();
		break;
	}

	const mat3d Fio = getRefConf()->getFinv();

	const double lt = (F * (Fio * curLayer->getBaseVector()[1])).norm();
	const double lz = (F * (Fio * curLayer->getBaseVector()[2])).norm();
	const double ore = curLayer->getTurnOverConstituents().at(0)->getReor();
	
	double n_angle = atan(tan(orig_angle) * pow(lt / lz, ore));
	curLayer->getTurnOverConstituents().at(0)->setAngle(n_angle);
	return n_angle;
}
bool Conf::recalcMassFrac(Layer* hLayer)
{
	try
	{

		/*for (size_t i = 0; i < m_layers.size(); i++)
		{*/
			Layer* refLayer = m_conf_ref->getCurrentLayer();
			double phi_ref = refLayer->getTurnOverConstituents().at(0)->getFraction();
			double phi_current = hLayer->getTurnOverConstituents().at(0)->getFraction();
			
			double eta = hLayer->getTurnOverRatio();

			double j_to_J = 1;

			switch (m_mode)
			{
			case MODE_POST:

				phi_current = getHomConf()->getCurrentLayer()->getTurnOverConstituents().at(0)->getFraction();
				j_to_J = m_conf_hom->getJ() / m_conf_ref->getJ();
				eta = getHomConf()->getCurrentLayer()->getTurnOverRatio();
				break;
			case MODE_PROCESS:
				//j_to_J = this->getJ() / m_conf_prv->getJ();
				j_to_J = this->getJ() / m_conf_ref->getJ();
				break;
			}


			
			double base = j_to_J * phi_current / phi_ref;

			if (m_mode == MODE_PROCESS)
			{
				phi_current = phi_ref;

				const double eps = std::numeric_limits<double>::epsilon();//0.000000000001;

				double dRdc = 1.0;
				double Rphi = 1.0;
				if (eta != 0)
				{
					while (true)
					{
						double ratio_sum = 1;
						double turn_over_res = 0;
						double no_turn_over_res = 0;

						for (size_t j = 1; j < refLayer->getTurnOverConstituents().size(); j++)
						{
							base = j_to_J * phi_current / phi_ref;
							double c_phi_ref = refLayer->getTurnOverConstituents()[j]->getFraction();

							ratio_sum += c_phi_ref / phi_ref * eta * pow(base, eta - 1.0);
							turn_over_res += c_phi_ref * pow(base, eta);
						}
						turn_over_res += j_to_J * (phi_current - 1);

						for (size_t j = 0; j < refLayer->getNoTurnOverConstituents().size(); j++)
						{
							no_turn_over_res += refLayer->getNoTurnOverConstituents()[j]->getFraction();
						}
						dRdc = j_to_J * ratio_sum;
						Rphi = turn_over_res + no_turn_over_res;
						phi_current -= Rphi / dRdc;

						if (abs(Rphi) <= sqrt(eps))
						{
							break;
						}
					}
				}
			}
			hLayer->getTurnOverConstituents()[0]->setMassFraction(phi_current);

			base = j_to_J * phi_current / phi_ref;
			//update turn over constituents mass fraction:
			for (size_t j = 1; j < hLayer->getTurnOverConstituents().size(); j++)
			{
				double phi = refLayer->getTurnOverConstituents()[j]->getFraction() / j_to_J * pow(base, eta);
				hLayer->getTurnOverConstituents()[j]->setMassFraction(phi);	
			}
			//update no turn over constituents mass fraction:
			for (size_t j = 0; j < hLayer->getNoTurnOverConstituents().size(); j++)
			{
				double phi = refLayer->getNoTurnOverConstituents()[j]->getFraction() / j_to_J;
				hLayer->getNoTurnOverConstituents()[j]->setMassFraction(phi);
			}
		return true;
	}
	catch (...)
	{
		std::string strStress712 = "D:/Error_in_recalc_mass_fraction.txt";
		std::ofstream file11(strStress712);

		return false;
	}
}
void Conf::updateMassFrac(Conf* pOther)
{
	for (std::size_t i = 0; i < m_layers.size(); i++)
	{

		// Update NoTurnOverConstituents
		for (std::size_t j = 0; j < getLayerByIndex(i)->getNoTurnOverConstituents().size(); j++)
		{
			Constituent* source = getLayerByIndex(i)->getNoTurnOverConstituents().at(j);
			Constituent* dest = pOther->getLayerByIndex(i)->getNoTurnOverConstituents().at(j);
			dest->setFraction(source->getFraction());
		}

		// Update TurnOverConstituents
		for (std::size_t j = 0; j < getLayerByIndex(i)->getTurnOverConstituents().size(); j++)
		{
			Constituent* source = getLayerByIndex(i)->getTurnOverConstituents().at(j);
			Constituent* dest = pOther->getLayerByIndex(i)->getTurnOverConstituents().at(j);
			dest->setFraction(dest->getFraction());
		}
	}
}
void Conf::updateConstituents(Conf* pOther)
{
}
void Conf::calcConstituentLevelConstantCauchyStress(mat3ds s, tens4dmm& stiffness_matrix)
{
	mat3ds	sfpro;
	vec3d	eigenvec[3];
	vec3d	Fxeigenvec[3];

	double	eigenval[3];

	sfpro.zero();
	m_U.eigen2(eigenval, eigenvec);

	sfpro(0, 0) = eigenvec[0] * (s * eigenvec[0]);
	sfpro(1, 1) = eigenvec[1] * (s * eigenvec[1]);
	sfpro(2, 2) = eigenvec[2] * (s * eigenvec[2]);
	sfpro(0, 1) = eigenvec[0] * (s * eigenvec[1]);
	sfpro(1, 2) = eigenvec[1] * (s * eigenvec[2]);
	sfpro(0, 2) = eigenvec[0] * (s * eigenvec[2]);

	Fxeigenvec[0] = m_F * eigenvec[0];
	Fxeigenvec[1] = m_F * eigenvec[1];
	Fxeigenvec[2] = m_F * eigenvec[2];

	
	for (int i = 0; i < 3; i++)
	{
		mat3ds ten1 = dyad(Fxeigenvec[i]);

		for (int j = 0; j < 3; j++)
		{

			double component = sfpro(i, j) / pow(eigenval[i], 3) / eigenval[j];

			mat3ds ten2 = dyads(Fxeigenvec[i], Fxeigenvec[j]);

			stiffness_matrix -= component * dyad1mm(ten2, ten1);

			for (int k = 0; k < 3; k++)
			{

				if (k == i) continue;

				mat3ds ten3 = dyads(Fxeigenvec[j], Fxeigenvec[k]);
				mat3ds ten4 = dyads(Fxeigenvec[k], Fxeigenvec[i]);

				component = sfpro(i, j) / eigenval[i] / eigenval[j] / eigenval[k] / (eigenval[i] + eigenval[k]);

				stiffness_matrix -= component * dyad1mm(ten3, ten4);
			}
		}
	}
}

double Conf::getCurrentLayerPenalty()
{
	setIrIo();
	const double inf = getCurrentLayer()->calcInfPn();
	const double flow1 = getCurrentLayer()->calcFlowPn();
	const double pn = 1 / (1.0 - getCurrentLayer()->getMechanoSensing());
	return pn * (1.0 + flow1 - inf);
}
tens4dmm Conf::getCurrentLayerPenalty(mat3ds& st, tens4dmm& stiffness)
{
	vec3d N[3]; getCurrentLayer()->getBaseVector(N);

	const mat3dd	I(1.0);
	const tens4ds	IxI = dyad1s(I);
	const tens4ds	IoI = dyad4s(I);
	const tens4dmm	Ixnrr = dyad1mm(I, dyad(m_F * (getRefConf()->getFinv() * N[0]))); // const
	const tens4dmm	Ixntt = dyad1mm(I, dyad(m_F * (getRefConf()->getFinv() * N[1]))); // const
	const tens4dmm	Ixsx = dyad1mm(I, st); // const
	const tens4dmm	IxIss = tens4dmm(IxI);							// IxI in tens4dmm form
	const tens4dmm	IoIss = tens4dmm(IoI);							// IoI in tens4dmm form

	//const double lr = (m_F * (getPrvConf()->getFinv() * N[0])).norm();						// lr -> 1 for F -> Fo
	//const double lt = (m_F * (getPrvConf()->getFinv() * N[1])).norm();
	const double lr = (m_F * (getRefConf()->getFinv() * N[0])).norm();						// lr -> 1 for F -> Fo
	const double lt = (m_F * (getRefConf()->getFinv() * N[1])).norm();
	const double rr = getCurrentLayer()->getRRatio();
	const double flow1 = getCurrentLayer()->calcFlowPn();
	const double flow2 = getCurrentLayer()->calcFlowDrvPn();
	const double inf = getCurrentLayer()->calcInfPn();
	//const double vv = this->getPrvConf()->getVolStress() / (1.0 - getCurrentLayer()->getMechanoSensing());
	const double vv = getVolStress() / (1.0 - getCurrentLayer()->getMechanoSensing());
	
	//return 1.0 / 3.0 * (2.0 * st.tr() * IoIss - 2.0 * Ixsx - ddot(IxIss, stiffness));
	return 1.0 / 3.0 * (2.0 * st.tr() * IoIss - 2.0 * Ixsx - ddot(IxIss, stiffness)) + vv * (1.0 + flow1 - inf) * (IxIss - 2.0 * IoIss) + vv * (flow2) * ((rr / lt) * Ixntt - ((rr - 1) / lr) * Ixnrr);
}

mat3ds Conf::calcStress()
{

	mat3ds s, S, Sx;

	mat3ds pn(0.0);
	mat3ds total_passive_turn(0.0);
	mat3ds total_passive_no_turn(0.0);
	mat3ds total_active_turn(0.0);
	mat3ds total_active_no_turn(0.0);
	double Jh_to_Jo = 1;
	static int x = 0;
	if (m_mode == MODE_PROCESS)
	{
		for (std::size_t i = 0; i < getCurrentLayer()->getNoTurnOverConstituents().size(); i++)
		{
			Constituent* constituent = getCurrentLayer()->getNoTurnOverConstituents()[i];
			constituent->init(this);
			double frc_org = getRefConf()->getCurrentLayer()->getNoTurnOverConstituents()[i]->getFraction();

			constituent->calcPassiveStress(this);
			total_passive_no_turn += frc_org * constituent->getPassiveCauchyStress();

			if (constituent->isActiveStress())
			{
				constituent->calcActiveStress(this);
				total_active_no_turn += frc_org * (constituent->getActiveCauchyStress(this));
			}
		}
		for (Constituent* constituent : getCurrentLayer()->getTurnOverConstituents())
		{
			constituent->init(this);
			double frc = constituent->getFraction();
//			if (frc != 0.0)
			{
				#pragma omp critical
				{
					constituent->calcPassiveStress(this);
				}
				if (frc != 0.0)
				{
					total_passive_turn += frc * (constituent->getPassiveCauchyStress());
					if (constituent->isActiveStress())
					{
						constituent->calcActiveStress(this);
						total_active_turn += frc * (constituent->getActiveCauchyStress(this));
					}
				}
			}
		}

		mat3d  ui(this->getUinv()); // converting mat3ds type to a mat3d
		Sx = total_passive_no_turn + total_active_no_turn + m_J * (ui * total_passive_turn * ui).sym() + m_J * (ui * total_active_turn * ui).sym();
		const double p = 1.0 / 3.0 / m_J * Sx.dotdot(m_C) - getVolStress() * getCurrentLayerPenalty();
		//return Sx;
		S = Sx - m_J * p * m_Ci;

	}
	else
	{
		//	case MODE_PRE/MODE_POST:
		for (std::size_t i = 0; i < getCurrentLayer()->getNoTurnOverConstituents().size(); i++)
		{
			Constituent* constituent = getCurrentLayer()->getNoTurnOverConstituents()[i];
			constituent->init(this);
			double frc = constituent->getFraction();
			if (m_mode == MODE_POST)
			{
				frc = getRefConf()->getCurrentLayer()->getNoTurnOverConstituents()[i]->getFraction();
			}
			total_passive_no_turn += (frc * (constituent->calcPassiveStress(this)));
			pn += m_Ci * constituent->calcScalarPenalty(this);

			if (constituent->isActiveStress() && m_mode == MODE_PRE) // probably defined as false
			{
				total_active_no_turn += (frc * (constituent->calcActiveStress(this)));
			}
		}

		for (std::size_t i = 0; i < getCurrentLayer()->getTurnOverConstituents().size(); i++)
		{
			Constituent* constituent = getCurrentLayer()->getTurnOverConstituents()[i];
	
			constituent->init(this);
			double frc = constituent->getFraction();
//			if (frc != 0.0)
			{
				total_passive_turn += (frc * (constituent->calcPassiveStress(this)));

				if (constituent->isActiveStress() && m_mode == MODE_PRE)
				{
					total_active_turn += (frc * (constituent->calcActiveStress(this)));
				}
			}
		}
		
		Sx = total_passive_no_turn + total_active_no_turn + total_passive_turn + total_active_turn;
		S = Sx + pn;

		m_volumetric_stress = 1.0 / 3.0 / m_J * (S).dotdot(m_C);
		//For calibration
		m_circ_stress	=	S.yy();
		m_axial_stress	=	S.zz();
	}

	s = 1.0 / m_J * ((m_F * ((S)* m_F.transpose()))).sym();

	return s;
}
tens4dmm Conf::calcTangent()
{
	tens4dmm s_tang(0.0);

	// Define identity tensor and some useful dyadic products of the identity tensor
	const mat3dd  I(1.0);
	const tens4ds IxI = dyad1s(I);
	const tens4ds IoI = dyad4s(I);
	const tens4dmm IxIss = tens4dmm(IxI);							// IxI in tens4dmm form
	const tens4dmm IoIss = tens4dmm(IoI);							// IoI in tens4dmm form

	if (m_mode == MODE_PROCESS)
	{
		//for (std::size_t layerIndex = 0; layerIndex < m_layers.size(); layerIndex++)
		//{
		//	Layer* layerPtr = (m_layers)[layerIndex];

			if (getCurrentLayer())
			{
				mat3d  ui(m_Ui);
				mat3ds Sx, sx, cnst_s;
				mat3ds tot_cnt_stress(0.0);
				mat3ds total_turn_over_cauchy_stress(0.0);
				mat3ds total_no_turn_over_cauchy_stress(0.0);
				mat3ds aS(0.0);

				tens4dmm	cass(0.0);
				tens4dmm	cess(0.0);
				tens4dmm	cpnss(0.0);
				tens4dmm	scxI2(0.0);

				

				for (std::size_t i = 0; i < getCurrentLayer()->getNoTurnOverConstituents().size(); i++)
				{
					Constituent* constituent = getCurrentLayer()->getNoTurnOverConstituents()[i];
					//double frc = constituent->getFraction();
					double frc_org = getRefConf()->getCurrentLayer()->getNoTurnOverConstituents()[i]->getFraction();
					mat3ds pS = constituent->getPassiveCauchyStress();
					total_no_turn_over_cauchy_stress += frc_org * pS;
					cess += frc_org * constituent->calcTangentPassiveComponent(this);
				}
				
				//	Defining the terms for collagen

				const auto& constituents = getCurrentLayer()->getTurnOverConstituents();
				double tot_non_col_mass_frc = 0;

				
				tens4dmm cfss(0.0);
				Constituent* constituent = constituents[0];
				double frc_col = constituent->getFraction();
				double frc_col_r = getRefConf()->getCurrentLayer()->getConstituents()[0]->getFraction();
////				
				mat3ds pS = m_J * (ui * constituent->getPassiveCauchyStress() * ui).sym();
				tens4dmm scxI = dyad1mm(1.0 / m_J * (m_F * (pS * m_F.transpose())).sym(), I);
				

				total_turn_over_cauchy_stress += frc_col * pS;
				tot_cnt_stress += frc_col * constituent->getPassiveCauchyStress();

				cpnss += constituent->calcTangetOrientationComponent(this);

				for (std::size_t i = 1; i < constituents.size(); i++)
				{
					Constituent* constituent = constituents[i];
					double frc = constituent->getFraction();
					double eta = getCurrentLayer()->getTurnOverRatio();
					double frc_r = getRefConf()->getCurrentLayer()->getConstituents()[i]->getFraction();
					double c_frc = frc_r * eta * pow(m_J / this->getRefConf()->getJ() * frc_col / frc_col_r, eta - 1.0);
					tot_non_col_mass_frc += c_frc;
				}
				for (std::size_t i = 1; i < constituents.size(); i++)
				{
					tens4dmm psxI(0.0);
					Constituent* constituent = constituents[i];

					double frc = constituent->getFraction();
			//		if (frc != 0)
					{
						tot_cnt_stress += frc * constituent->getPassiveCauchyStress();
						
							mat3ds pS = m_J * (ui * constituent->getPassiveCauchyStress() * ui).sym();
						if (frc != 0)
						{
							total_turn_over_cauchy_stress += frc * pS;

							psxI += dyad1mm(1.0 / m_J * (m_F * (pS * m_F.transpose())).sym(), I);
						}
						//////
						//double eta = getCurrentLayer()->getTurnOverRatio();
						////double frc_r = this->getPrvConf()->getLayerByIndex(layerIndex)->getConstituents()[i]->getFraction();
						//double frc_r = getRefConf()->getCurrentLayer()->getConstituents()[i]->getFraction();
						//double c_frc = frc_r * eta * pow(m_J / this->getRefConf()->getJ() * frc_col / frc_col_r, eta - 1.0);
						//tot_non_col_mass_frc += c_frc;

						if (constituent->isActiveStress())
						{
							tot_cnt_stress += frc * constituent->getActiveCauchyStress(this);
							aS = m_J * (ui * constituent->getActiveCauchyStress(this) * ui).sym();
							psxI += dyad1mm(1.0 / m_J * (m_F * (aS * m_F.transpose())).sym(), I); // const
							total_turn_over_cauchy_stress += frc * aS;

							cass += frc * constituent->calcTangentActiveComponent(this);
						}

						double eta = getCurrentLayer()->getTurnOverRatio();
						double frc_r = getRefConf()->getCurrentLayer()->getConstituents()[i]->getFraction();
						double c_frc = frc_r * eta * pow(m_J / this->getRefConf()->getJ() * frc_col / frc_col_r, eta - 1.0);

						double dphiRm = c_frc / (tot_non_col_mass_frc + frc_col_r); // const
						cfss += dphiRm * (psxI);
					}
				}
				

				calcConstituentLevelConstantCauchyStress(tot_cnt_stress, cfss);

				double dphiRm = tot_non_col_mass_frc / (tot_non_col_mass_frc + frc_col_r); // const
				double dphiRc = 1 - dphiRm; // const
				
				//cfss += dphiRm * (psxI) + dphiRc * scxI;
				cfss += dphiRc * scxI;
	
				

				mat3ds total_s = total_no_turn_over_cauchy_stress + total_turn_over_cauchy_stress;
				sx = 1.0 / m_J * (m_F * (total_s * m_F.transpose())).sym();
				tens4dmm css = cess + cfss + cass + cpnss;
				s_tang += css + getCurrentLayerPenalty(sx, css);
			}
		//}
		
	}
	else
		//case MODE_PRE / MODE_POST
	{
		for (std::size_t i = 0; i < getCurrentLayer()->getNoTurnOverConstituents().size(); i++)
		{
			Constituent* constituent = getCurrentLayer()->getNoTurnOverConstituents()[i];
			double frc = constituent->getFraction();

			s_tang += frc * (constituent->calcTangentPassiveComponent(this));

			s_tang += constituent->calcTangentPenalty(this);

			if (constituent->isActiveStress() && m_mode == MODE_PRE) // Probably false
			{
				s_tang += frc * (constituent->calcTangentActiveComponent(this)) / m_J;
			}
		}

		for (Constituent* constituent : m_currentLayer->getTurnOverConstituents())
		{

			double frc = constituent->getFraction();
			s_tang += frc * (constituent->calcTangentPassiveComponent(this) / m_J);


			if (constituent->isActiveStress() && m_mode == MODE_PRE)
			{
				s_tang += frc * (constituent->calcTangentActiveComponent(this)) / m_J;
			}
		}
	}

	return (s_tang);
}
Conf& Conf::operator=(const Conf& other)
{
	m_F = other.m_F;
	m_C = other.m_C;
	m_J = other.m_J;
	m_current_t = other.m_current_t;
	m_Ci = other.m_Ci;
	m_Fi = other.m_Fi;
	m_U = other.m_U;
	m_Ui = other.m_Ui;
	m_R = other.m_R;
	m_volumetric_stress = other.m_volumetric_stress;
	m_conf_ref = nullptr;
	m_conf_hom = nullptr;
	m_conf_prv = nullptr;
	m_currentLayer = other.m_currentLayer;
	m_partial_t = other.m_partial_t;
	m_mode = other.m_mode;

	// Copying layers data
	//m_layers.clear();
	//
	//for (const Layer* layerPtr : other.m_layers)
	//{
	//	m_layers.push_back(new Layer(*layerPtr));
	//}

	return *this;
}