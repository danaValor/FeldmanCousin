//simpson intergration
#include<iostream>
#include<cmath>
using namespace std;



double mu(double x, double y, double E1,double E2, double L1, double L2)
{
	auto ff = [&](double E) {
		double r = (L2 - L1 - 1 / (2 * 1.27 * y) * E * (sin(2 * 1.27 * y * L2 / E) - sin(2 * 1.27 * y * L1 / E))) / 2 / (L2 - L1);
		return r;
	};



	auto Simpson = [&](double a, double b, int n) {//a is lower bound, b is upper bound
	
		double h = (b - a) / n;
		double x = a, sum0 = 0, sum = 0;
		for (int i = 0; i < n; i++)
		{
			if ((i == 0£© || £¨i = n - 1£©) sum0 = ff(x) * h / 3;
			else if (i % 2 = 1) sum0 = ff(x) * 2 * h / 3;
			else sum0 = ff(x) * 4 * h / 3;
			sum += sum0;
			x += a;
		}

		return sum;
	};


}








}