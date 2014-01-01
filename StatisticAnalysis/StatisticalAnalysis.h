#pragma once

#include <string>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <limits>
#include <iomanip>
#include <fstream>



	class StatisticalAnalysis
	{
	public:
		
		double FDistributionInverse(double probability, int m, int n);

		double FDistributionSearch(double probability, int m, int n, int step, double start, double end);

		double FDistribution(double x, int freedom1, int freedom2);

		double BetaCF(double a, double b, double x);

		double BetaFunction(double m, double n);

		double BetaIncomplete(double a, double b, double x);

		double GammLn(double n);
	};
