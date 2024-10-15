#ifndef _INSULT_H
#define _INSULT_H
#define _USE_MATH_DEFINES
#include <math.h>
#include <algorithm>

#include "Layer.h"
#include "Constituent.h"

#endif // !INSULT_H
class Layer;
/*
The isult class reperesents a pertubation, which is defined by its type, value, spatial params and the layer.
Since some of the insult are specific to a given constituent, dynamic cast was used.
*/
enum InsultType
{
	//General Insult
	FlowInsult,
	//
	InflammationRatioInsult,
	InflammationAmplitudeInsult,
	//
	MechanoSensingInsult,
	GeneralProductionRemovalInsult, //eta
	//Constituent Insult
	ElasticFiberInsult,
	//
	SmContractilityInsult,			//Tmax
	//
	CollagenCrossLinkingInsult,		//cc
	CollagenMechanoRegulationInsult //Gc
};
class Insult
{
public:
	Insult(InsultType insType, double insult, Layer* pLocation);
	void setAxialParam(double slope,double width, double apex);
	void setCircParam(double slope, double width, double apex, double asym);
	void setInsultPeriod(double start, double end);
	void applyInsult(Layer* pUpdate);
private:
	//properties
	InsultType m_insType;
	double m_insult;
	Layer* layer;

	//factors
	double m_timeFactor;
	double m_axialFactor;
	double m_circumFactor;
};
class InsManager 
{
public:
	InsManager() {}
	~InsManager() {}
	void registerInsult(const Insult& insult)
	{
		m_registeredInsults.push_back(insult);
	}
	void applyInsults(Layer* u)
	{

		for (auto& insult : m_registeredInsults)
		{
			insult.applyInsult(u);
		}
	}
	void UpdateInsultsPeriod(double current_time, double partial_time)
	{
		for (auto& insult : m_registeredInsults)
		{
			insult.setInsultPeriod(current_time, partial_time);
		}
	}
private:
	std::vector<Insult> m_registeredInsults;

};

