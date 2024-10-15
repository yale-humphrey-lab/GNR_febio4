#ifndef _CONSTITUENT_H
#define _CONSTITUENT_H

#include "FECore/mat3d.h"						// To get current time
#include "Conf.h"

class Conf;
class Layer;
/*
This is a pure virtual class to define new constituents. The computation of stress is based on a given configuration.
One can extend the level of abstruction by defining the constituent laws in a configuration file.
In addition a fibre class can be defined and utlise for various constituents.
*/
class Constituent
{
public:
	Constituent();
	Constituent(double mass, double modulus, double G, bool to, double angle, double reorient);
	Constituent(double mass, double modulus, double G, bool to, double reorient);
	virtual Constituent* clone() const = 0;
	virtual ~Constituent();
	// Energy and Stress Methods
	void setParameters(double mass, double modulus, double G, bool to, double reorient);
	void setParameters(double mass, double modulus, double G, bool to, double angle, double reorient);
	virtual void init(Conf* pCnf) = 0;
	
	virtual double calcScalarPenalty(Conf* pCnf) = 0;

	virtual bool isPenalty() = 0;
	virtual bool isActiveStress() = 0;
	virtual bool isTurnOver() = 0;
	
	virtual void calcEnergy(Conf* pCnf) = 0;
	virtual mat3ds calcActiveStress(Conf* pCnf) =  0;
	virtual mat3ds calcPassiveStress(Conf* pCnf) = 0;

	
	virtual tens4dmm calcTangentActiveComponent(Conf* pCnf) = 0;
	virtual tens4dmm calcTangentPenalty(Conf* pCnf) = 0;
	virtual tens4dmm calcTangetOrientationComponent(Conf* pCnf) = 0;
	virtual tens4dmm calcTangentPassiveComponent(Conf* pCnf) = 0;


	//Get
	mat3ds getActiveCauchyStress(Conf* pCnf, bool decay = false) { return m_act_cauchy_stress; }
	mat3ds getPassiveCauchyStress() { return m_psv_cauchy_stress; }

	void recalcAngle(Conf* pCnf);

	double getFraction();
	void setFraction(const double frc) { m_mass_frac = frc; }
	double getReor() { return m_reorientation; }
	
	//Set
	void setMassFraction(double frac);
	void setAngle(double alpha);
	double getAngle() { return m_angle; }

	void setModulus(double m) { m_modulus = m; }
	const double getModulus() const { return m_modulus; }
	void setStretch(double G) { m_G = G; }
	const double getStretch() const { return m_G; }

protected:
	void calcCauchyStress(Conf* pCnf) {};

	
	//Mass fraction
	
	
	double m_angle;		// alphao: Initial diagonal collagen orientation (rad)
	double m_modulus;	// Modulus a/Shear
	double m_mass_frac;

	double m_G;
	double m_lm;
	//double m_energy;

	bool   m_turn_over;
	double m_reorientation;

	mat3ds m_psv_cauchy_stress;
	mat3ds m_act_cauchy_stress;
	
};

#endif // _CONSTITUENT_H
