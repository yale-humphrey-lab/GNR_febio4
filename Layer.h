#ifndef _LAYER_H
#define _LAYER_H
#include <vector>
#include "Constituent.h"

#include <mutex>

class Constituent;

/*
	This class defines a layer with different constituents by implementing a polymorphic approach.
	Due to legacy reasons, the constituents are stored in two separate data structures, to diffrentiate between turn-over and non turn-over behaviour.
	Using the isTurnOver method, one can easily merge the two vectors into a single one.
	Please note that both the RRatio and Flow Penalties assume a cylindrical tube.
*/
class Layer
{
public:
	Layer();
	Layer(double inner, double lo, double to_ratio = 1.0);
	Layer(const Layer& other);
	~Layer();

	void load(double inner, double length, double to_ratio = 1.0);
	//Constituents methods
	void addConstituent(Constituent* constituent);

	const std::vector<Constituent*>& getConstituents() const;
	const std::vector<Constituent*>& getTurnOverConstituents() const;
	const std::vector<Constituent*>& getNoTurnOverConstituents() const;

	const double getTurnOverRatio() { return turnover_ratio; }
	void		 setTurnOverRatio(double r) { turnover_ratio = r; }
	

	//Geometrical dimensions
	void setPosition(vec3d& p);
	void init(double inner, double length);
	const double getInnerRadius() { return rIo; }
	const double getAxialLength() { return lo; }
	const double getRRatio() { return r_ratio; }
	void setRRatio(double r) { r_ratio = r / rIo; }
	void setIrIo(double lt, double lr) { rIrIo = r_ratio * lt - (r_ratio - 1) * lr; }
	void		 getBaseVector(vec3d* other);
	const vec3d* getBaseVector() const { return N; }
	void setBaseVec(vec3d x, vec3d y, vec3d z) { N[0] = x; N[1] = y; N[2] = z; }
	const double getZ() { return material_point.z; }
	void setAzimuth(double azm) { m_azi = azm; }
	const double getAzimuth() const { return m_azi; }
	const double getAxialLength() const { return lo; }
	//characteristics
	const double getMechanoSensing() const { return m_mechanosensing; }
	void setMechanoSensing(double sn) { m_mechanosensing = sn; }
	void setTortuosity(double amnt, double prd) { m_tortuosity_amount = amnt; m_tortuosity_period = prd; }
	//Penalties 
	const double getFlowPn() const{ return flow_pn; }
	void setFlowPn(double f) { flow_pn = f; }
	const double calcFlowPn() { return flow_pn * (pow(rIrIo, -3) - 1.0); }
	const double calcFlowDrvPn() { return -3 * flow_pn * pow(rIrIo, -4); }
	//
	const double getInfPn() { return inflammatory_amp_pn * inflammatory_ratio_pn; }
	void setInfAmpPn(double inf) { inflammatory_amp_pn = inf; }
	void setInfRatioPn(double inf) { inflammatory_ratio_pn = inf; }
	const double calcInfPn() { return inflammatory_amp_pn * inflammatory_ratio_pn; }

private:

	double rIo;         // Initial inner radius
	double lo;          // Initial axial length
	double r_ratio;			// Radial vector size
	double rIrIo;
	vec3d N[3];			// Base vector
	vec3d material_point;
	double m_azi;			// Azimuth wrt axis -Y

	//penalties
	double flow_pn;
	double inflammatory_amp_pn;
	double inflammatory_ratio_pn;
	
	double radial_def_pn;

	//characteristics
	double m_tortuosity_amount;			// Tortuosity amount
	double m_tortuosity_period;			// Tortuosity period
	double m_mechanosensing;			// Delta
	//constituents members
	double turnover_ratio;// production_removal_turnover: eta;
	std::vector<Constituent*> m_turn_over_constituents;
	std::vector<Constituent*> m_no_turn_over_constituents;
	mutable std::mutex LayerMutex_;
public:
	Layer& operator=(const Layer& other);
};

#endif