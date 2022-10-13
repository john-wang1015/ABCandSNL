#pragma once

#define MAX_DBL std::numeric_limits<double>::max()

inline double Sqr(double x)
{
	return x * x;
}

inline double DistanceSquared(double X1, double Y1, double X2, double Y2)
{
	double dx = X1 - X2;
	double dy = Y1 - Y2;
	return Sqr(dx) + Sqr(dy);
}

