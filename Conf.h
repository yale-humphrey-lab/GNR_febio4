#ifndef _CONF_H
#define _CONF_H

#include "FECore/FEAnalysis.h"					// To get end time
#include "FECore/FEModel.h"						// To get current time
#include "FECore/mat3d.h"
#include "FECore/tens4d.h"
#include <vector>
#include <memory>
#include <mutex>
#include "Constituent.h"
#include "Layer.h"

/*
	This class models a specific configuration and calculates the composition and stress magnatudes of the constituents.
	The class includes a vector of layers, which can be define both in the radial and axial direction.
	The first can be utilised to characterise several layers of the vessel e.g., media & adventitia in the artery.
	The later can be utilised to characterise end-to-end graft anastomosis.
	Originally, the Conf class held a pointer to a vector of Layers, but due to multithreading limitations of FEBio, these data structures were externalised.
*/

class Layer;
class Constituent;

enum Mode 
{
	MODE_PRE = 0,
	MODE_PROCESS = 1,
	MODE_POST = 2
};
class Conf
{
public:
	Conf();
	Conf(const Conf& other);
	~Conf();
	void setMode();
	const Mode getMode();
	void init(mat3d& F, mat3ds& C, double J, double t);
	int setIrIo();
	int setPosition(vec3d& v, bool tor = false);
	///////////////////////////////////////
	void addLayer(Layer* input_layer);
	void   setCurrentLayer(Layer* input_layer);
	Layer* getCurrentLayer();
	Layer* getCurrentLayer(vec3d& v);
	Layer* getLayerByIndex(std::size_t layerIndex);
	//////////////////////////////////////
	double getVolStress() { return m_volumetric_stress; }
	double getStressComp(int num) { if (num == 3) { return m_axial_stress; } else { return m_circ_stress; }}

	const mat3d& getR() { return m_R; }
	const mat3d& getF() { return m_F; }
	const mat3d& getFinv() { return m_Fi; }
	const double getJ() { return m_J; };
	const mat3ds& getU() { return m_U; }
	const mat3ds& getUinv() { return m_Ui; }

	Conf* getRefConf() { return m_conf_ref; }
	Conf* getPrvConf() { return m_conf_prv; }
	Conf* getHomConf() { return m_conf_hom; }
	/////////////////////////////////////////////
	void setRefConf(Conf* pOther);
	void setHomConf(Conf* pOther);
	void setPreviousConf(Conf* pOther);

	double	recalcFibersAngle(Layer* curLayer);
	bool	recalcMassFrac(Layer* hLayer);
	void	updateMassFrac(Conf* pOther);
	void	updateConstituents(Conf* pOther);
	void	calcConstituentLevelConstantCauchyStress(mat3ds s, tens4dmm& stiffness_matrix);
	double	getCurrentLayerPenalty();
	tens4dmm	getCurrentLayerPenalty(mat3ds& st, tens4dmm& stiffness);
	/////////////////////////////////////////////
	mat3ds	 calcStress();
	tens4dmm calcTangent();
	/////////////////////////////////////////////
	Conf& operator=(const Conf& other);
private:
	//Refrence to other Conf instances
// Constructor

	std::vector<Layer*> m_layers;
	Conf* m_conf_ref;
	Conf* m_conf_prv; // if accounting more than the ref and hom configs.
	Conf* m_conf_hom;

	Layer* m_currentLayer;

	//Mode
	Mode m_mode;

	double m_volumetric_stress;

	double m_current_t;
	double m_partial_t;

	mat3d m_F;
	mat3d m_Fi;
	double m_J;

	mat3d  m_R;
	mat3ds m_U;
	mat3ds m_Ui;

	mat3ds m_C;
	mat3ds m_Ci;

	double m_circ_stress;
	double m_axial_stress;
};
#endif // !CONF_H
