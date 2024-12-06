#include "FEMbeCmm.h"
#include "FEBioMech/FEElasticMaterial.h"
#include "FECore/FEDomain.h"
#include "FECore/FEAnalysis.h"					// to get end time
#include "FECore/FEModel.h"						// to get current time
#include "FECore/log.h"							// to print to log file and/or screen
#include <iostream>								// to use cin.get()
#include <sstream>
#include <signal.h>
#define _USE_MATH_DEFINES						// to introduce pi constant (1/2)
#include <math.h>								// to introduce pi constant (2/2)
#include <limits>
///
#include <ctime>
#include <fstream>
#include <cstdlib>

// define the material parameters
BEGIN_FECORE_CLASS(FEMbeCmm, FEElasticMaterial)
/*
partial time - to define the scope of the G&R
end time - additional cycle with elevated pressure - usually unnecessary
*/
ADD_PARAMETER(endtime, FE_RANGE_GREATER(0.0), "endtime");
ADD_PARAMETER(partialtime, FE_RANGE_GREATER(0.0), "partialtime");
ADD_PARAMETER(m_layer_id, "layerid");
END_FECORE_CLASS();

FEMaterialPointData* GRMaterialPoint::Copy()
{
	GRMaterialPoint* pt = new GRMaterialPoint(*this);
	if (m_pNext) pt->m_pNext = m_pNext->Copy();
	return pt;
}

void GRMaterialPoint::Init()
{
//	ref_config = nullptr;
//	last_config = nullptr;
//	current_config = nullptr;

//	FEMaterialPointData::Init();
//	InitLayers();
	RegisterInsults();
}

void GRMaterialPoint::Serialize(DumpStream& ar)
{
	FEMaterialPointData::Serialize(ar);
	//ar& m_Jo& m_svo& m_smo& m_sco& m_Fio& m_Jh& m_Fih& m_phic& m_Iemax& m_stress_inv_h;
}
///////////
void GRMaterialPoint::InitLayers()
{
	/*
	This is a naive parser to collect the data of each layer and constituents.
	It is required that the first turn over constituent within each layer would be collagen.
	*/
	FEMaterialPoint& pt = *(this->ExtractData<FEMaterialPoint>());

	std::ifstream inputFile("D:/noG/materials.txt"); // Change to local directory
	//std::ifstream inputFile("materials.txt");

	if (!inputFile.is_open())
	{
		std::string Strrsm = "D:/Cannot open materials file";
		std::ofstream rf(Strrsm);
		return;
	}
	int numLayers = 0;
	
	std::string line;

	while (std::getline(inputFile, line))
	{
		if (line.find("Layer:") != std::string::npos)
		{
			numLayers++;
		}
	}

	inputFile.clear();  // Reset the file stream state
	inputFile.seekg(0); // Rewind to the beginning of the file

	// Initialize an array of Layer pointers with the calculated size
	m_orig_layers.reserve(numLayers);
	m_prv_layers.reserve(numLayers);
	m_current_layers.reserve(numLayers);

	while (std::getline(inputFile, line))
	{
		if (line.find("Layer:") != std::string::npos)
		{
			int layer_indx = 0;
			std::istringstream lineStream(line);

			string type;
			double rIo, lo;
			double deltab = 0.0;
			double KsKib = 0.0;
			double etab = 1.0;
			double imper = 0.0;
			double hwaves = 2.0;
			lineStream >> type >> rIo >> lo >> deltab >> KsKib >> etab >> imper >> hwaves;
			Layer* layer = new Layer(rIo, lo, etab);
			layer->setMechanoSensing(deltab);
			layer->setFlowPn(KsKib);
			layer->setTortuosity(imper, hwaves);

			//layer->setPosition(dPoint);

			std::getline(inputFile, line);
			std::istringstream eStream(line);

			double phieo, Get, Gez, mub, Jdep, bulkLM;
			eStream >> type >> phieo >> Get >> Gez >> mub >> Jdep >> bulkLM;
			Elastin* els = new Elastin(phieo, Get, Gez, mub, Jdep, bulkLM);
			layer->addConstituent(els);

			std::getline(inputFile, line);
			std::istringstream cStream(line);
			double phico, Gcb, ccb, dc, betat, betaz, alphao, aexp;
			cStream >> type >> phico >> Gcb >> ccb >> dc >> betat >> betaz >> alphao >> aexp;

			Collagen* collagen = new Collagen(phico, Gcb, ccb, dc, betat, betaz, alphao, aexp);
			layer->addConstituent(collagen);

			std::getline(inputFile, line);
			std::istringstream sStream(line);
			double phimo, Gm, cm, dm, Tmaxb, CB, lamM, lam0;
			if (sStream >> type >> phimo >> Gm >> cm >> dm >> Tmaxb >> CB >> lamM >> lam0)
			{
				if (type == "Sm:")
				{
					Sm* sm = new Sm(phimo, Gm, cm, dm, Tmaxb, CB, lamM, lam0);
					layer->addConstituent(sm);
				}
			}

			std::getline(inputFile, line);
			std::istringstream gStream(line);

			double phigo, gGrad, gub, gJdep, gbulkLM;
			if (gStream >> type >> phigo >> gGrad >> gub >> gJdep >> gbulkLM)
			{
				if (type == "GaGs:")
				{
					Gags* gag = new Gags(phigo, gGrad, gub, gJdep, gbulkLM);
					layer->addConstituent(gag);
				}
			}

			m_orig_layers.push_back(layer);
		}
	}

	inputFile.close();

	ref_config = new Conf();

	for (const Layer* layerPtr : m_orig_layers)
	{
		m_prv_layers.push_back(new Layer(*layerPtr));
		m_current_layers.push_back(new Layer(*layerPtr));
	}

	last_config = new Conf();
	current_config = new Conf();

}

void GRMaterialPoint::RegisterInsults()
{
	std::ifstream inputFile("D:/noG/Insult.txt"); //Change to local folder
	//std::ifstream inputFile("Insult.txt"); //Change to local folder

	if (!inputFile.is_open())
	{
		std::string Strrsm = "D:/noG/Cannot open materials file";
		std::ofstream rf(Strrsm);
		return;
	}

	std::string line;

	while (std::getline(inputFile, line))
	{
		double insult;
		char delimiter;
		int layer_index;
		std::string type;
		std::istringstream lineStream(line);
		lineStream >> type >> insult >> layer_index >> delimiter;

		InsultType insultType;
		if (type == "CollagenCrossLinkingInsult") {
			insultType = InsultType::CollagenCrossLinkingInsult;
		}
		else if (type == "ElasticFiberInsult") {
			insultType = InsultType::ElasticFiberInsult;
		}
		else if (type == "SmContractilityInsult") {
			insultType = InsultType::SmContractilityInsult;
		}
		else if (type == "CollagenMechanoRegulationInsult") {
			insultType = InsultType::CollagenMechanoRegulationInsult;
		}
		else if (type == "MechanoSensingInsult") {
			insultType = InsultType::MechanoSensingInsult;
		}
		else if (type == "GeneralProductionRemovalInsult") {
			insultType = InsultType::GeneralProductionRemovalInsult;
		}
		else if (type == "FlowInsult") {
			insultType = InsultType::FlowInsult;
		}
		else if (type == "InflammationAmplitudeInsult") {
			insultType = InsultType::InflammationAmplitudeInsult;
		}
		else if (type == "InflammationRatioInsult") {
			insultType = InsultType::InflammationRatioInsult;
		}
		else {
			// Unknown insult type
			std::string strStressr742 = "D:/Unknown type" + type;
			std::ofstream rfil4114(strStressr742);
			return;
		}
		int prm_cnt = 0;
		Insult ins(insultType, insult, m_orig_layers.at(layer_index));
		// Parse and initialize axial parameters
		double nuz, zod, zo;
		if (lineStream >> nuz >> zod >> zo)
		{
			ins.setAxialParam(nuz, zod, zo);
		}
		// Parse and initialize circular parameters
		double nut, tod, to, asym;
		if (lineStream >> delimiter >> nut >> tod >> to >> asym)
		{
			ins.setCircParam(nut, tod, to, asym);
		}

		m_Insults.registerInsult(ins);
	}

	inputFile.close();
}
///////////


FEMaterialPointData* FEMbeCmm::CreateMaterialPointData()
{
	return new GRMaterialPoint(new FEElasticMaterialPoint);
}

FEMbeCmm::FEMbeCmm(FEModel* pfem) : FEElasticMaterial(pfem)
{
	m_secant_tangent = true;
}

bool FEMbeCmm::Init()
{
	if (!FEElasticMaterial::Init()) return false;

	FEDomainList& domList = GetDomainList();

	for (int i = 0; i < domList.Domains(); ++i)
	{

		FEDomain& dom = *domList.GetDomain(i);


		for (int j = 0; j < dom.Elements(); ++j)

		{
			FEElement& el = dom.ElementRef(j);
			for (int n = 0; n < el.GaussPoints(); ++n)
			{

				FEMaterialPoint& mp = *el.GetMaterialPoint(n);
				GRMaterialPoint* gr = mp.ExtractData<GRMaterialPoint>();

				if (gr)
				{
					gr->ref_config = nullptr;
					gr->last_config = nullptr;
					gr->current_config = nullptr;
						// TODO: initialize layers here instead of GRMaterialPoint::Init
					std::ifstream inputFile("D:/noG/materials.txt"); // Change to local directory
					//std::ifstream inputFile("materials.txt");

					if (!inputFile.is_open())
					{
						std::string Strrsm = "D:/Cannot open materials file";
						std::ofstream rf(Strrsm);
						return false;
					}
					int numLayers = 0;
					std::string line;

					while (std::getline(inputFile, line))
					{
						if (line.find("Layer:") != std::string::npos)
						{
							numLayers++;
						}
					}

					inputFile.clear();  // Reset the file stream state
					inputFile.seekg(0); // Rewind to the beginning of the file

					// Initialize an array of Layer pointers with the calculated size
					gr->m_orig_layers.reserve(numLayers);
					gr->m_prv_layers.reserve(numLayers);
					gr->m_current_layers.reserve(numLayers);

					while (std::getline(inputFile, line))
					{
						if (line.find("Layer:") != std::string::npos)
						{
							int layer_indx = 0;
							std::istringstream lineStream(line);

							string type;
							double rIo, lo;
							double deltab = 0.0;
							double KsKib = 0.0;
							double etab = 1.0;
							double imper = 0.0;
							double hwaves = 2.0;
							lineStream >> type >> rIo >> lo >> deltab >> KsKib >> etab >> imper >> hwaves;
							Layer* layer = new Layer(rIo, lo, etab);
							layer->setMechanoSensing(deltab);
							layer->setFlowPn(KsKib);
							layer->setTortuosity(imper, hwaves);
							layer->setPosition(mp.m_r0);

							std::getline(inputFile, line);
							std::istringstream eStream(line);

							double phieo, Get, Gez, mub, Jdep, bulkLM;
							eStream >> type >> phieo >> Get >> Gez >> mub >> Jdep >> bulkLM;
							Elastin* els = new Elastin(phieo, Get, Gez, mub, Jdep, bulkLM);
							layer->addConstituent(els);

							std::getline(inputFile, line);
							std::istringstream cStream(line);
							double phico, Gcb, ccb, dc, betat, betaz, alphao, aexp;
							cStream >> type >> phico >> Gcb >> ccb >> dc >> betat >> betaz >> alphao >> aexp;

							Collagen* collagen = new Collagen(phico, Gcb, ccb, dc, betat, betaz, alphao, aexp);
							layer->addConstituent(collagen);

							std::getline(inputFile, line);
							std::istringstream sStream(line);
							double phimo, Gm, cm, dm, Tmaxb, CB, lamM, lam0;
							if (sStream >> type >> phimo >> Gm >> cm >> dm >> Tmaxb >> CB >> lamM >> lam0)
							{
								if (type == "Sm:")
								{
									Sm* sm = new Sm(phimo, Gm, cm, dm, Tmaxb, CB, lamM, lam0);
									layer->addConstituent(sm);
								}
							}

							std::getline(inputFile, line);
							std::istringstream gStream(line);

							double phigo, gGrad, gub, gJdep, gbulkLM;
							if (gStream >> type >> phigo >> gGrad >> gub >> gJdep >> gbulkLM)
							{
								if (type == "GaGs:")
								{
									Gags* gag = new Gags(phigo, gGrad, gub, gJdep, gbulkLM);
									layer->addConstituent(gag);
								}
							}

							gr->m_orig_layers.push_back(layer);
						}
					}

					inputFile.close();

					gr->ref_config = new Conf();

					for (const Layer* layerPtr : gr->m_orig_layers)
					{
						gr->m_prv_layers.push_back(new Layer(*layerPtr));
						gr->m_current_layers.push_back(new Layer(*layerPtr));
					}

					gr->last_config = new Conf();
					gr->current_config = new Conf();


					double layer_index = m_layer_id(mp);

					gr->ref_config->setCurrentLayer(gr->m_orig_layers.at(layer_index));
					gr->last_config->setCurrentLayer(gr->m_prv_layers.at(layer_index));
					gr->current_config->setCurrentLayer(gr->m_current_layers.at(layer_index));
				}

			}

		}

	}
////
}

void FEMbeCmm::StressTangent(FEMaterialPoint& mp, mat3ds& stress, tens4dmm& tangent)
{
	// The FEMaterialPoint classes are stored in a linked list. The specific material
	// point data needed by this function can be accessed using the ExtractData member.
	// In this case, we want to get FEElasticMaterialPoint data since it stores the deformation
	// information that is needed to evaluate the stress.
	FEElasticMaterialPoint& et = *mp.ExtractData<FEElasticMaterialPoint>();
	GRMaterialPoint& pt = *mp.ExtractData<GRMaterialPoint>();

	// We'll need the deformation gradient and its determinant in this function.
	// Note that we don't take the determinant of F directly (using mat3d::det)
	// but instead use the m_J member variable of FEElasticMaterialPoint.

	const double eps = std::numeric_limits<double>::epsilon();

	// Get current and end times
	const double t = GetFEModel()->GetTime().currentTime;
	const double sgr = min(t, partialtime);

	mat3d& F = et.m_F;
	double J = et.m_J;
	mat3ds C = et.RightCauchyGreen();

	//In case there are more than one layer - define different mp according to spatial distrib.
	double layer_index = m_layer_id(mp);

	//After pre-stretch
	if (t > 1 + eps)
	{
		//Some insult are time dependent.
		pt.m_Insults.UpdateInsultsPeriod(t, partialtime);
		//Apply insults by layer
		pt.m_Insults.applyInsults(pt.m_current_layers.at(layer_index));
	}
	//Pre-stretch step
	if (t <= 1 + eps)
	{
		pt.ref_config->init(F, C, J, t);
		pt.current_config->init(F, C, J, t);

		stress = pt.current_config->calcStress();
		tangent = pt.current_config->calcTangent();

		*pt.last_config = *pt.current_config;
		pt.last_config->setCurrentLayer(pt.m_prv_layers.at(layer_index));

		return;
	}
	//Hom step
	else if (t <= partialtime + eps)
	{
		pt.current_config->init(F, C, J, t);

		pt.current_config->setRefConf(pt.ref_config);
		pt.current_config->setPreviousConf(pt.last_config);

		pt.current_config->recalcFibersAngle(pt.m_current_layers.at(layer_index));

		pt.current_config->recalcMassFrac(pt.m_current_layers.at(layer_index));

		pt.current_config->setCurrentLayer(pt.m_current_layers.at(layer_index));

		stress = pt.current_config->calcStress();
		tangent = pt.current_config->calcTangent();


		for (std::size_t i = 0; i < pt.m_current_layers.size(); i++)
		{
			*pt.m_prv_layers[i] = *pt.m_current_layers[i];
		}
		*pt.last_config = *pt.current_config;
		pt.last_config->setCurrentLayer(pt.m_prv_layers.at(layer_index));


		if (layer_index == 0)
		{
			et.m_a.x = pt.current_config->getCurrentLayer()->getNoTurnOverConstituents().at(0)->getFraction();
			et.m_a.y = pt.current_config->getCurrentLayer()->getTurnOverConstituents().at(0)->getFraction();
			et.m_a.z = pt.current_config->getCurrentLayer()->getTurnOverConstituents().at(1)->getFraction();
			double ggs = 1 - (et.m_a.x + et.m_a.y + et.m_a.z);
			feLog("FMF_media elastin %.4f\t collagen %.4f\t sm % .4f\t GaGs %.4f\n", et.m_a.x, et.m_a.y, et.m_a.z, ggs);
		}
		else
		{
			et.m_v.x = pt.current_config->getCurrentLayer()->getNoTurnOverConstituents().at(0)->getFraction();
			et.m_v.y = pt.current_config->getCurrentLayer()->getTurnOverConstituents().at(0)->getFraction();
			et.m_v.z = pt.current_config->getCurrentLayer()->getTurnOverConstituents().at(1)->getFraction();
			double ggs = 1 - (et.m_v.x + et.m_v.y + et.m_v.z);
			feLog("FMF_adv elastin %.4f\t collagen %.4f\t sm % .4f\t GaGs %.4f\n", et.m_v.x, et.m_v.y, et.m_v.z, ggs);
		}

		return;
	}

	//Post hom - G&R is inactive - usually unnecessary
	else {

		if (layer_index == 0)
		{
			et.m_a.x = pt.current_config->getCurrentLayer()->getNoTurnOverConstituents().at(0)->getFraction();
			et.m_a.y = pt.current_config->getCurrentLayer()->getTurnOverConstituents().at(0)->getFraction();
			et.m_a.z = pt.current_config->getCurrentLayer()->getTurnOverConstituents().at(1)->getFraction();
			double ggs = 1 - (et.m_a.x + et.m_a.y + et.m_a.z);
			feLog("FMF_media elastin %.4f\t collagen %.4f\t sm % .4f\t GaGs %.4f\n", et.m_a.x, et.m_a.y, et.m_a.z, ggs);
		}
		else
		{
			et.m_v.x = pt.current_config->getCurrentLayer()->getNoTurnOverConstituents().at(0)->getFraction();
			et.m_v.y = pt.current_config->getCurrentLayer()->getTurnOverConstituents().at(0)->getFraction();
			et.m_v.z = pt.current_config->getCurrentLayer()->getTurnOverConstituents().at(1)->getFraction();
			double ggs = 1 - (et.m_v.x + et.m_v.y + et.m_v.z);
			feLog("FMF_adv elastin %.4f\t collagen %.4f\t sm % .4f\t GaGs %.4f\n", et.m_v.x, et.m_v.y, et.m_v.z, ggs);
		}


		{
			pt.last_config->setRefConf(pt.ref_config);

			pt.current_config->init(F, C, J, t);
			pt.current_config->setRefConf(pt.ref_config);
			pt.current_config->setHomConf(pt.last_config);

			pt.last_config->recalcFibersAngle(pt.m_current_layers.at(layer_index));

			pt.current_config->recalcMassFrac(pt.m_current_layers.at(layer_index));

			pt.current_config->setCurrentLayer(pt.m_current_layers.at(layer_index));

			stress = pt.current_config->calcStress();
			tangent = pt.current_config->calcTangent();
		}
	}
}
