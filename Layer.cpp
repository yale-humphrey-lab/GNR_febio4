#include "Layer.h"
#include "Insult.h"
#include "Constituent.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>
#include <fstream>
#include <sstream>

Layer::Layer()
{
	m_azi = 0;
	m_tortuosity_amount = 0;
	m_mechanosensing = 0;
	//Penalties:
	flow_pn = 0.35;
	inflammatory_amp_pn = 0;
	inflammatory_ratio_pn = 0;
	m_turn_over_constituents.clear();
	m_no_turn_over_constituents.clear();
}
Layer::Layer(double in, double lng, double to_ratio)
{
	rIo = in;
	lo = lng;

	m_azi = 0;
	m_tortuosity_amount = 0;
	m_mechanosensing	= 0;
	
	turnover_ratio		= to_ratio;

	flow_pn = 0.35;
	inflammatory_amp_pn = 0;
	inflammatory_ratio_pn = 0;
	m_turn_over_constituents.clear();
	m_no_turn_over_constituents.clear();
}
Layer::Layer(const Layer& other)
{
	m_turn_over_constituents.clear();
	m_no_turn_over_constituents.clear();
	
	for (const Constituent* ptr : other.getNoTurnOverConstituents())
	{
		if (ptr)
		{
			m_no_turn_over_constituents.push_back(ptr->clone());

		}
		else
		{
			std::string strStress7712 = "D:/addConstituent no turn is null.txt";
			std::ofstream file41777772211(strStress7712);
		}
	}

	for (const Constituent* ptr : other.getTurnOverConstituents())
	{
		if (ptr)
		{
			m_turn_over_constituents.push_back(ptr->clone());

		}
		else
		{
			std::string strStress7712 = "D:/addConstituent turn over is null.txt";
			std::ofstream file41777772211(strStress7712);
		}
	}

	rIo = other.rIo; lo = other.lo; r_ratio = other.r_ratio;
	m_tortuosity_amount = other.m_tortuosity_amount;
	m_tortuosity_period = other.m_tortuosity_period;
	inflammatory_amp_pn = other.inflammatory_amp_pn;
	inflammatory_ratio_pn = other.inflammatory_ratio_pn;
	radial_def_pn = other.radial_def_pn;
	m_mechanosensing = other.m_mechanosensing;		
	turnover_ratio = other.turnover_ratio;
	m_azi = other.m_azi;
	flow_pn = other.flow_pn; 
	rIrIo = other.rIrIo;
	N[0] = other.N[0]; N[1] = other.N[1]; N[2] = other.N[2];
	material_point = other.material_point;

}
Layer::~Layer()
{
	for (Constituent* constituent_to : m_turn_over_constituents)
		delete constituent_to;
	for (Constituent* constituent_n_to : m_turn_over_constituents)
		delete constituent_n_to;
	m_turn_over_constituents.clear();
	m_no_turn_over_constituents.clear();
}
void Layer::load(double inner, double length, double to_ratio)
{
	rIo = inner;
	lo = length;

	m_azi = 0;
	m_tortuosity_amount = 0;
	m_mechanosensing = 0;
	inflammatory_amp_pn = 0;
	inflammatory_ratio_pn = 0;

	turnover_ratio = to_ratio;

	flow_pn = 0;

	m_turn_over_constituents.clear();
	m_no_turn_over_constituents.clear();
}
void Layer::init(double inner, double length)
{
	rIo = inner; lo = length;
}
////////////////////////////////////////////////////////////////////////////////
void Layer::addConstituent(Constituent* constituent)
{
	if (!constituent)
	{
		std::string strStress7712 = "D:/Constituent is null cannot add.txt";
		std::ofstream file41777772211(strStress7712);
		return;
	}
	
	if (constituent->isTurnOver()) 
	{
		m_turn_over_constituents.push_back(constituent);
	}
	else 
	{
		m_no_turn_over_constituents.push_back(constituent);	
	}
}
const std::vector<Constituent*>& Layer::getConstituents() const
{
	vector<Constituent*> AllConstitPtr = m_turn_over_constituents; // Copy the first vector to the result
	AllConstitPtr.insert(AllConstitPtr.end(), m_no_turn_over_constituents.begin(), m_no_turn_over_constituents.end()); // Append the second vector
	return AllConstitPtr;
}
const std::vector<Constituent*>& Layer::getNoTurnOverConstituents() const
{
	return m_no_turn_over_constituents;
}
const std::vector<Constituent*>& Layer::getTurnOverConstituents() const 
{
	return m_turn_over_constituents;
}
////////////////////////////////////////////////////////////////////////////////
void Layer::setPosition(vec3d& p)
{
	material_point = p;

	double R = (2 * lo) / M_PI;


	double phi = atan2(p.y, p.z);
	double theta = atan2(p.x, (sqrt(pow(p.y, 2.) + pow(p.z, 2.)) - R));

	if (true)
	{
											
		const vec3d E_r = { sin(theta), cos(theta) * sin(phi), cos(theta) * cos(phi)};
		const vec3d E_theta = { cos(theta), -sin(theta) * sin(phi), -sin(theta) * cos(phi)};
		const vec3d E_phi = { 0, -cos(phi), sin(phi)};

		setBaseVec(E_r, E_theta, E_phi);

		//double coef = 5;
		double coef = 5;
		double loc = 0.5; //0.2-0.5
		//double width = 4.5; 4- 5.5 if injury 0.2 and not 0.3 then 2.0-7.5
		double width = 7.5;

		double tu = 1/width * M_PI/4;

		m_azi = exp(pow(-fabs((phi - loc * M_PI/2) / tu), coef));
		//std::string strStressr742 = "D:/m_azi" + std::to_string(m_azi);
		//std::ofstream rfil4114(strStressr742);
	}
	else 
	{
		const vec3d Xcl = { 0.0, m_tortuosity_amount / 100.0 * rIo * sin(m_tortuosity_period * M_PI * p.z / lo), p.z };		// Center line
		vec3d NX = { p.x - Xcl.x, p.y - Xcl.y, p.z - Xcl.z };								// Radial vector
		double ro = sqrt(NX * NX);
		NX /= ro; //Er
		setRRatio(ro);

		vec3d N0, N1, N2;
		N2 = { 0.0, m_tortuosity_amount / 100.0 * rIo * m_tortuosity_period * M_PI / lo * cos(m_tortuosity_period * M_PI * p.z / lo), 1.0 };
		N2 /= sqrt(N2 * N2);  // Ez
		N1 = { -NX.y, NX.x, NX.z };		// Etheta																			// Circumferential
		N0 = N2 ^ N1; // Ez X Etheta = Er = Nx

		setBaseVec(N0, N1, N2);

		m_azi = acos(-NX.y);
	}		
}
void Layer::getBaseVector(vec3d* other)
{
	for (int i = 0; i < 3; ++i) 
	{
		other[i] = N[i];
	}	
}
Layer& Layer::operator=(const Layer& other)
{
	// Clear existing data
	m_turn_over_constituents.clear();
	m_no_turn_over_constituents.clear();

	// Copy constituents
	for (Constituent* ptr : other.getConstituents())
	{
		this->addConstituent(ptr->clone());
	}

	// Copy other members
	rIo = other.rIo;
	lo = other.lo;
	r_ratio = other.r_ratio;
	m_tortuosity_amount = other.m_tortuosity_amount;
	m_tortuosity_period = other.m_tortuosity_period;
	inflammatory_amp_pn = other.inflammatory_amp_pn;
	inflammatory_ratio_pn = other.inflammatory_ratio_pn;
	radial_def_pn = other.radial_def_pn;
	m_mechanosensing = other.m_mechanosensing;
	turnover_ratio = other.turnover_ratio;
	m_azi = other.m_azi;
	flow_pn = other.flow_pn;
	rIrIo = other.rIrIo;
	N[0] = other.N[0];
	N[1] = other.N[1];
	N[2] = other.N[2];
	material_point = other.material_point;
	
	return *this;
}
//////////////////////////////////////////////////////////////////////////////