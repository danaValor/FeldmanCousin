
#define _CRT_SECURE_NO_WARNINGS
#include<TString.h>
#include<iostream>
using namespace std;

int main()
{

		int const Nn = 100, bin = 5;
		int nGenerated[bin] = { 95, 100, 130, 125, 110 };
		int b[bin] = { 100 };
		double mu[bin] = { 100, 110, 120, 130, 140 };
		double chi = 0;
		for (int i = 0; i < 5; i++) {
			double chi0 = (nGenerated[i] - b[i] - mu[i])*(nGenerated[i] - b[i] - mu[i]) / (b[i] + mu[i]);
			chi = chi + chi0;
		}


		cout << chi << endl;
		getchar();
		return 0;


	TString a = "Hellow, ROOT";
	printf("%s\n", (char const *)a);
	return 0;

}
