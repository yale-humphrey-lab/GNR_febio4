//#pragma once
//=============================================================================
// This plugin example illustrates how to create a new material. 
// It requires FEBio 2.5 (or up)
//
// Author Steve Maas
// Copyright (c) 2015 - 2016
// All rights reserved
//
//=============================================================================

//-----------------------------------------------------------------------------
// We need to include this file since our new material class will inherit from
// FEElasticMaterial which is defined in this include files.
#include "FEBioMech/FEElasticMaterial.h"
#include <iostream>								// to use cin.get()

#include <list>
#include "Sm.h"
#include "Elastin.h"
#include "Collagen.h"
#include "GaGs.h"
#include "Layer.h"
#include "Conf.h"
#include "Insult.h"

class GRMaterialPoint : public FEMaterialPointData
{
public:
	GRMaterialPoint(FEMaterialPointData* pt) : FEMaterialPointData(pt) {};

	FEMaterialPointData* Copy() override;

	void Init() override;
	void Serialize(DumpStream& ar) override;

	void InitLayers();

	void RegisterInsults();

	~GRMaterialPoint()
	{
		for (Layer* lyr : m_orig_layers)
		{
			for (Constituent* cnt : lyr->getNoTurnOverConstituents())
			{
				delete cnt;
				cnt = nullptr;
			}
			for (Constituent* cnt : lyr->getTurnOverConstituents())
			{
				delete cnt;
				cnt = nullptr;
			}
			m_orig_layers.clear();
		}
		for (Layer* lyr : m_prv_layers)
		{
			for (Constituent* cnt : lyr->getNoTurnOverConstituents())
			{
				delete cnt;
				cnt = nullptr;
			}
			for (Constituent* cnt : lyr->getTurnOverConstituents())
			{
				delete cnt;
				cnt = nullptr;
			}
			m_prv_layers.clear();
		}
		for (Layer* lyr : m_current_layers)
		{
			for (Constituent* cnt : lyr->getNoTurnOverConstituents())
			{
				delete cnt;
				cnt = nullptr;
			}
			for (Constituent* cnt : lyr->getTurnOverConstituents())
			{
				delete cnt;
				cnt = nullptr;
			}
			m_current_layers.clear();
		}
		if (ref_config)
		{
			delete ref_config;
			ref_config = nullptr;
		}
	}

public:

	Conf* ref_config;
	Conf* last_config;
	Conf* current_config;

	InsManager m_Insults; // This is a factory class to define various pertubations

	//The Layers are defined here due to limitation of FEBio

	std::vector<Layer*> m_orig_layers;
	std::vector<Layer*> m_prv_layers;
	std::vector<Layer*> m_current_layers;
};

//-----------------------------------------------------------------------------
// This material class implements a neo-Hookean constitutive model. 
// Since it is a (compressible, coupled) hyper-elatic material, it needs to inherit
// from FEElasticMaterial. 
class FEMbeCmm : public FEElasticMaterial
{
public:
	FEMbeCmm(FEModel* pfem);
	bool Init();
	FEMaterialPointData* CreateMaterialPointData() override;

public:
	// The constructor is called when an instance of this class is created.
	// All classes registered by the framework must take the FEModel* as the only
	// parameter in the constructor, even if the class does not need it (which most often
	// will be the case). For material classes, the FEModel parameter is passed to the 
	// base class in the initialization list.

	// setting m_secant_tangent = true so FESolidMaterial uses SecantTangent
	// (allows minor symmetry only tangents) instead of Tangent (minor and major symmetries)
	// 	 { return m_secant_tangent; }
	bool m_secant_tangent = true;   //!< flag for using secant tangent
	double endtime;     // Simulation end time
	double partialtime; // G&R end time
	FEParamDouble m_layer_id; //MP for different layers

	DECLARE_FECORE_CLASS();

public:
	// function to perform material evaluation. calculates stress and tangent to avoid code duplication
	void StressTangent(FEMaterialPoint& mp, mat3ds& stress, tens4dmm& tangent);

	// This function calculates the spatial (i.e. Cauchy or true) stress.
	// It takes one parameter, the FEMaterialPoint and returns a mat3ds object
	// which is a symmetric second-order tensor.
	virtual mat3ds Stress(FEMaterialPoint& pt) override {
		mat3ds stress;
		tens4dmm tangent;
		StressTangent(pt, stress, tangent);
		return stress;
	}

	// This function calculates the spatial elasticity tangent tensor. 
	// It takes one parameter, the FEMaterialPoint and retursn a tens4ds object
	// which is a fourth-order tensor with major and minor symmetries.
	virtual tens4ds Tangent(FEMaterialPoint& pt) override {
		tens4ds tangent;
		return tangent;
	};

	// minor symmetries only
	virtual tens4dmm SecantTangent(FEMaterialPoint& pt) {
		mat3ds stress;
		tens4dmm tangent;
		StressTangent(pt, stress, tangent);
		return tangent;
	}

	bool UseSecantTangent() override { return m_secant_tangent; }

};
