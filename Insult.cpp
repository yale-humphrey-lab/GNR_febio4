#include "Insult.h"
#include "Sm.h"
#include "Elastin.h"
#include "Collagen.h"
#include <fstream>
#include <math.h>

Insult::Insult(InsultType insType, double insult, Layer* pLocation) : 
	m_insType(insType), m_insult(insult), layer(pLocation), m_timeFactor(1.0), m_axialFactor(1.0), m_circumFactor(1.0)
{
}
void Insult::setAxialParam(double slope, double width, double apex)
{
	//m_axialFactor = exp(-pow(fabs((layer->getZ() - apex) / (layer->getAxialLength() / 2.0 / width)), slope));
	m_axialFactor = layer->getAzimuth();
}
void Insult::setCircParam(double slope, double width, double apex, double asym)
{
	m_circumFactor = exp(-asym * pow(abs((layer->getAzimuth() - apex * M_PI) / (M_PI / width)), slope));
}
void Insult::setInsultPeriod(double start, double end)
{
	m_timeFactor = (min(start, end) - 1.0) / (end - 1.0);
}
 void Insult::applyInsult(Layer* pUpdate)
 {
	 int index		= -1;
	 double calc	= 0.0;
	 double orig	= 0.0;
	 const double insult_factor = m_axialFactor * m_circumFactor * m_timeFactor;
	 
	 std::string Strrsm;

	 switch (m_insType)
	 {
	 case	FlowInsult:
		 orig = layer->getFlowPn();
		 pUpdate->setFlowPn(orig + (m_insult - orig) * m_axialFactor);
		 break;
	 case InflammationAmplitudeInsult:
		 pUpdate->setInfAmpPn(m_insult); //* insult_factor);
		 break;
	 case	InflammationRatioInsult:
		 pUpdate->setInfRatioPn(m_insult * m_axialFactor);
		 break;
	 case	MechanoSensingInsult: //delta
		 orig = layer->getMechanoSensing();
		 pUpdate->setMechanoSensing(orig + ((m_insult) * insult_factor));
		 break;
	 case	GeneralProductionRemovalInsult: //eta
		 orig = layer->getTurnOverRatio();
		 pUpdate->setTurnOverRatio(orig + (m_insult - orig) * m_axialFactor);
		 break;
	 case	ElasticFiberInsult:
		 index = -1;

		 for (int i = 0; i < layer->getNoTurnOverConstituents().size(); i++) 
		 {
			 Constituent* c = layer->getNoTurnOverConstituents()[i];
			 if (Elastin * Ptr = dynamic_cast<Elastin*>(c)) 
			 {
				 orig = Ptr->getModulus();
				 index = i;
				 break;
			 }
		 }
		 if (index != -1) 
		 {	 

			 Constituent* c = pUpdate->getNoTurnOverConstituents()[index];
			 if (Elastin * Ptr = dynamic_cast<Elastin*>(c)) 
			 {
				// Ptr->setModulus(orig * (1 - m_insult * insult_factor));
				 Ptr->setModulus(orig * (1 - m_insult * insult_factor));
			 }
		 }
		 break;
	 case	SmContractilityInsult:			 //Tmax
		index = -1;

		 for (int i = 0; i < layer->getTurnOverConstituents().size(); i++)
		 {
			 Constituent* c = layer->getTurnOverConstituents()[i];
			 if (Sm * Ptr = dynamic_cast<Sm*>(c))
			 {
				 orig = Ptr->getActiveParam();
				 index = i;
				 break;
			 }
		 }
		 if (index != -1)
		 {
			 Constituent* c = pUpdate->getTurnOverConstituents()[index];
			 if (Sm * Ptr = dynamic_cast<Sm*>(c))
			 {
				 Ptr->setActiveParam(orig * (1 - m_insult * m_axialFactor * m_timeFactor));
			 }
		 }
		 break;
	 case	CollagenCrossLinkingInsult:		 //cc
		index = -1;

		 for (int i = 0; i < layer->getTurnOverConstituents().size(); i++)
		 {
			 Constituent* c = layer->getTurnOverConstituents()[i];
			 if (Collagen * Ptr = dynamic_cast<Collagen*>(c))
			 {
				 orig = Ptr->getModulus();
				 index = i;
				 break;
			 }
		 }
		 if (index != -1)
		 {
			 Constituent* c = pUpdate->getTurnOverConstituents()[index];
			 if (Collagen * Ptr = dynamic_cast<Collagen*>(c))
			 {
				 Ptr->setModulus(orig * (1 - m_insult * insult_factor));
			 }
		 }
		 break;
	 case	CollagenMechanoRegulationInsult: //Gc
		index = -1;

		 for (int i = 0; i < layer->getTurnOverConstituents().size(); i++)
		 {
			 Constituent* c = layer->getTurnOverConstituents()[i];
			 if (Collagen * Ptr = dynamic_cast<Collagen*>(c))
			 {
				 orig = Ptr->getStretch();
				 index = i;
				 break;
			 }
		 }
		 if (index != -1)
		 {
			 Constituent* c = pUpdate->getTurnOverConstituents()[index];
			 if (Collagen * Ptr = dynamic_cast<Collagen*>(c))
			 {
				 Ptr->setStretch(orig * (1 - m_insult * insult_factor));
			 }
		 }
		 break;
	 }
 }