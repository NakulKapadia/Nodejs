#include "statisticalanalysis.h"
#include <iostream>
#include <string>
#include <cmath>
#include <limits>
#include <iomanip>
#include <fstream>
//Node Code
#include <node.h>
#include <cstdlib>
#include <ctime>
using namespace v8;

// This function returns a JavaScript number that is either 0 or 1.
Handle<Value> FINV(const Arguments& args) {
	// At the top of every function that uses anything about v8, include a
	// definition like this. It ensures that any v8 handles you create in that
	// function are properly cleaned up. If you see memory rising in your
	// application, chances are that a scope isn't properly cleaned up.
	HandleScope scope;

	// Check that there are enough arguments. If we access an index that doesn't
	// exist, it'll be Undefined().
	//if (args.Length() < 1) {
	//    // No argument was passed. Throw an exception to alert the user to
	//    // incorrect usage. Alternatively, we could just use 0.
	//    return ThrowException(
	//        Exception::TypeError(String::New("First argument must be a number"))
	//    );
	//}

	// Cast a value to a specific type. See
	// http://izs.me/v8-docs/classv8_1_1Value.html for available To*() functions
	// and type checking functions. When converting to integer, make sure the
	// POD type you use is big enough!
	Local<Number> v1 = Local<Number>::Cast(args[0]);
	Local<Number> v2 = Local<Number>::Cast(args[1]);
	Local<Number> v3 = Local<Number>::Cast(args[2]);
	double probability = ((double)v1->Value()) ;
	int m = ((int)v2->Value()) ;
	int n = ((int)v3->Value()) ;

	StatisticalAnalysis* statsana = new StatisticalAnalysis();
	double res = statsana->FDistributionInverse(probability,m,n);
	/*if (aa3 < 0) {
		return ThrowException(Exception::TypeError(String::New(
			"Fibonacci sequence number must be positive")));
	}*/
	return scope.Close(Local<Value>::New(Number::New(res)));
}


Handle<Value> FDIST(const Arguments& args) {
	// At the top of every function that uses anything about v8, include a
	// definition like this. It ensures that any v8 handles you create in that
	// function are properly cleaned up. If you see memory rising in your
	// application, chances are that a scope isn't properly cleaned up.
	HandleScope scope;

	// Check that there are enough arguments. If we access an index that doesn't
	// exist, it'll be Undefined().
	//if (args.Length() < 1) {
	//    // No argument was passed. Throw an exception to alert the user to
	//    // incorrect usage. Alternatively, we could just use 0.
	//    return ThrowException(
	//        Exception::TypeError(String::New("First argument must be a number"))
	//    );
	//}

	// Cast a value to a specific type. See
	// http://izs.me/v8-docs/classv8_1_1Value.html for available To*() functions
	// and type checking functions. When converting to integer, make sure the
	// POD type you use is big enough!
	Local<Number> v1 = Local<Number>::Cast(args[0]);
	Local<Number> v2 = Local<Number>::Cast(args[1]);
	Local<Number> v3 = Local<Number>::Cast(args[2]);
	double probability = ((double)v1->Value()) ;
	int m = ((int)v2->Value()) ;
	int n = ((int)v3->Value()) ;

	StatisticalAnalysis* statsana = new StatisticalAnalysis();
	double res = statsana->FDistribution(probability,m,n);
	/*if (aa3 < 0) {
		return ThrowException(Exception::TypeError(String::New(
			"Fibonacci sequence number must be positive")));
	}*/
	return scope.Close(Local<Value>::New(Number::New(res)));
}

void RegisterModule(Handle<Object> target) {
	// target is the module object you see when require()ing the .node file.
	target->Set(String::NewSymbol("FINV"),
		FunctionTemplate::New(FINV)->GetFunction());
	target->Set(String::NewSymbol("FDIST"),
		FunctionTemplate::New(FDIST)->GetFunction());
}

NODE_MODULE(StatisticalAnalysis, RegisterModule);


double StatisticalAnalysis::FDistributionInverse(double probability, int m, int n)
{
	if (probability == 0.0)
	{
		return std::numeric_limits<double>::quiet_NaN(); 
	}
	if (probability == 1.0)
	{
		return 0.0;
	}
	if ((probability < 0.0) || (probability > 1.0))
	{
		throw "ExceptionStatisticalAnalysesInvalidProbabilityValue";
	}
	int step = 0;
	return this->FDistributionSearch(probability, m, n, step, 0.0, 10000.0);
}

double StatisticalAnalysis::FDistributionSearch(double probability, int m, int n, int step, double start, double end)
{
	step++;
	double x = (start + end) / 2.0;
	double num2 = this->FDistribution(x, m, n);
	if (step > 30)
	{
		return x;
	}
	if (num2 <= probability)
	{
		return this->FDistributionSearch(probability, m, n, step, start, x);
	}
	return this->FDistributionSearch(probability, m, n, step, x, end);
}

double StatisticalAnalysis::FDistribution(double x, int freedom1, int freedom2)
{
	if (x < 0.0)
	{
		throw "ExceptionStatisticalAnalysesInvalidTValue";
	}
	if (freedom1 <= 0)
	{
		throw "ExceptionStatisticalAnalysesInvalidDegreeOfFreedom";
	}
	if (freedom2 <= 0)
	{
		throw "ExceptionStatisticalAnalysesInvalidDegreeOfFreedom";
	}
	if (x == 0.0)
	{
		return 1.0;
	}
	if (x == std::numeric_limits<double>::quiet_NaN())
	{
		return 0.0;
	}
	return this->BetaIncomplete((static_cast<double>(freedom2)) / 2.0, (static_cast<double>(freedom1)) / 2.0, (static_cast<double>(freedom2)) / (freedom2 + (freedom1 * x)));
}

double StatisticalAnalysis::BetaCF(double a, double b, double x)
{
	int num = 100;
	double num2 = 3E-07;
	double num3 = 1E-30;
	double num11 = a + b;
	double num13 = a + 1.0;
	double num12 = a - 1.0;
	double num7 = 1.0;
	double num8 = 1.0 - ((num11 * x) / num13);
	if (abs(num8) < num3)
	{
		num8 = num3;
	}
	num8 = 1.0 / num8;
	double num10 = num8;
	int num4 = 1;
	while (num4 <= num)
	{
		int num5 = 2 * num4;
		double num6 = ((num4 * (b - num4)) * x) / ((num12 + num5) * (a + num5));
		num8 = 1.0 + (num6 * num8);
		if (abs(num8) < num3)
		{
			num8 = num3;
		}
		num7 = 1.0 + (num6 / num7);
		if (abs(num7) < num3)
		{
			num7 = num3;
		}
		num8 = 1.0 / num8;
		num10 *= num8 * num7;
		num6 = ((-(a + num4) * (num11 + num4)) * x) / ((a + num5) * (num13 + num5));
		num8 = 1.0 + (num6 * num8);
		if (abs(num8) < num3)
		{
			num8 = num3;
		}
		num7 = 1.0 + (num6 / num7);
		if (abs(num7) < num3)
		{
			num7 = num3;
		}
		num8 = 1.0 / num8;
		double num9 = num8 * num7;
		num10 *= num9;
		if (abs(static_cast<double>(num9 - 1.0)) < num2)
		{
			break;
		}
		num4++;
	}
	if (num4 > num)
	{
		throw "ExceptionStatisticalAnalysesIncompleteBetaFunction";
	}
	return num10;
}

double StatisticalAnalysis::BetaFunction(double m, double n)
{
	return exp((this->GammLn(m) + this->GammLn(n)) - this->GammLn(m + n));
}

double StatisticalAnalysis::BetaIncomplete(double a, double b, double x)
{
	double num;
	if ((x < 0.0) || (x > 1.0))
	{
		throw"ExceptionStatisticalAnalysesInvalidInputParameter";
	}
	if ((x == 0.0) || (x == 1.0))
	{
		num = 0.0;
	}
	else
	{
		num = exp((((this->GammLn(a + b) - this->GammLn(a)) - this->GammLn(b)) + (a * log(x))) + (b * log(1.0 - x)));
	}
	if (x < ((a + 1.0) / ((a + b) + 2.0)))
	{
		return ((num * this->BetaCF(a, b, x)) / a);
	}
	return (1.0 - ((num * this->BetaCF(b, a, 1.0 - x)) / b));
}

double StatisticalAnalysis::GammLn(double n)
{
	double num;
	double numArray[6] = {76.180091729471457, -86.505320329416776, 24.014098240830911, -1.231739572450155, 0.001208650973866179, -5.395239384953E-06};
	if (n < 0.0)
	{
		throw "ExceptionStatisticalAnalysesGammaBetaNegativeParameters";
	}
	double num2 = num = n;
	double d = num + 5.5;
	d -= (num + 0.5) * log(d);
	double num4 = 1.0000000001900149;
	for (int i = 0; i <= 5; i++)
	{
		num4 += numArray[i] / ++num2;
	}
	return (-d + log((2.5066282746310007 * num4) / num));
}


