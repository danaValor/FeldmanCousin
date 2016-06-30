
#define _CRT_SECURE_NO_WARNINGS

#define R__EXTERN       R__DllImport extern

#include<TString.h>
#include<TRandom3.h>
#include<TMinuit.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<algorithm>
#include <math.h>
#include <stdio.h>
//#include<TFitterMinuit.h>
//#include<TMinuit.h>
//#include<TROOT.h>

using namespace std;

TRandom3 rnd(0);
int const bin = 5;
int const background = 100;
int const maxNn = 10000;
int Nn = 1000;
double Ebins[] = { 10, 20, 30, 40, 50, 60 };
double L1 = 0.6, L2 = 1;
static int nExperiment[bin];
static int *nexp = &nExperiment[0];
double initx = 0.0001;
double xlower = -4.;
int nx = 200, ny = 200;
int const loopx = 10, loopy = 10;

double Simpson(int n, double E1, double E2, double y) {//a is lower bound, b is upper bound

	double h = (E2 - E1) / n;
	double e = E1, sum0 = 0, sum = 0;
	for (int i = 0; i < n; i++) {
		double r = (L2 - L1 - 1 / (2 * 1.27 * y) * e * (sin(2 * 1.27 * y * L2 / e) -
			sin(2 * 1.27 * y * L1 / e))) / 2;
		if (i == 0 || i == n - 1) sum0 = r * h / 3;
		else {
			if (i % 2 == 1) sum0 = r * 2 * h / 3;
			else sum0 = r * 4 * h / 3;
		}
		sum += sum0;
		e += h;
	};


	return sum;
};


double Mu(double x, double y, int binNumber)
{
	double E2 = Ebins[binNumber + 1];
	double E1 = Ebins[binNumber];


	int const N = 3000;
	double muNumber = x * Simpson(N, E1, E2, y) * 10000 / (L2 - L1) / (E2 - E1);// n = 1000, which is the number of integration bins

	return muNumber;
}



int Ngenerated(double muIn)
{

	double sigma = sqrt(muIn + background);
	double mean = muIn + background;
	double nn = rnd.Gaus(mean, sigma);
	int ng = int(round(nn));

	return ng;

}




double Chi(double x, double y) //binNumber = 0.. 4
{
	double chiInit = 0;
	for (int i = 0; i < bin; i++) {
		chiInit += (nexp[i] - background - Mu(x, y, i)) *
			(nexp[i] - background - Mu(x, y, i)) / (background + Mu(x, y, i));

	}


	return chiInit;

}




//int nPar = 2;

void MinuitFunction(int& nDim, double* gout, double& result, double par[], int flg)
{
	//double a0 = 1.5;
	result = Chi(par[0], par[1]);

}


double FindMinimum(double xI, double yI)
{
	//TMinuit minimizer(2);

//	minimizer.mninit(5, 6, 7);
//	minimizer.SetPrintLevel(-1);
	int flag = 0;

//	minimizer.SetFCN(MinuitFunction);
//	minimizer.DefineParameter(0, "xx", xI, 0.01, 0, 1);
//	minimizer.DefineParameter(1, "yy", yI, 1, 0, 1000);
//	minimizer.mnexcm("SIMPLEX", 0, 0, flag);

//	minimizer.SetErrorDef(1);

//	minimizer.SetMaxIterations(500);
//	int err = minimizer.Migrad();
//	if (err != 0) {
//		return err;
//	}

//	double xx1, mu1err;
//	double yy1, mu2err;

//	minimizer.GetParameter(0, xx1, mu1err);
//	minimizer.GetParameter(1, yy1, mu2err);



	//after finding the global minimum with simplex, use migrad to find the local minimum.
	TMinuit minimizer1(2);
	minimizer1.mninit(5, 6, 7);
	minimizer1.SetPrintLevel(-1);
	minimizer1.SetFCN(MinuitFunction);
	minimizer1.DefineParameter(0, "xx2", xI, 0.01, 0, 1);
	minimizer1.DefineParameter(1, "yy2", yI, 1, 0, 1000);

	minimizer1.SetErrorDef(1);

	minimizer1.SetMaxIterations(1000);
	int err1 = minimizer1.Migrad();
	if (err1 != 0) {
		return 10000;
	}

	double xx21, mu1err1;
	double yy21, mu2err1;

	minimizer1.GetParameter(0, xx21, mu1err1);
	minimizer1.GetParameter(1, yy21, mu2err1);

	double minimum = Chi(xx21, yy21);

	return minimum;

}



double ChiCritical(double x, double y) {
	double chiCritical[loopx][loopy] = { 0 };
	double deltaChi[maxNn] = { 0 };

	double mu[bin] = { 0 };


	int cut = int(ceil(0.9 * Nn));


	for (int p = 0; p < Nn; p++) {
		for (int i = 0; i < bin; i++) {
			mu[i] = Mu(x, y, i);

			nExperiment[i] = Ngenerated(mu[i]);

		}
		double chiMin = FindMinimum(x, y);
		deltaChi[p] = Chi(x, y) - chiMin;
		if (deltaChi[p] < 0) {
			p = p - 1;
			cout << "find minimum of chi-square again\n";
		}

	}

	sort(deltaChi, deltaChi + Nn);
	double chiCri = deltaChi[cut];
	return chiCri;

}



void WriteIn(double x,double y)
{
	char str[100];
	sprintf(str, "data/x_%.4f_y_%.4f", x,y);
	ofstream f1(str);
	double chix1;
	chix1 = ChiCritical(x,y);
	f1 << chix1 << endl;
	f1.close();
}






int main(int argc, char *argv[])
{

	int i, j, n;
	
//	printf("%s", argv[0]);
	sscanf(argv[1], "%d", &i);
	sscanf(argv[2], "%d", &j);
	sscanf(argv[3], "%d", &n);
//	printf("%d %d %d", i, j, n);

	Nn = n;

	double x0 = initx, y0 = 1;

	double logx0 = log10(x0);
	double logy0 = log10(y0);
	double hx = (0 - logx0) / (loopx - 1);
	double hy = (3 - logy0) / (loopy - 1);
	double chiCritical[loopx][loopy] = { 0 };
	

	logx0 += i * hx;
	logy0 += j * hy;
	x0 = exp(log(10)*logx0);
	y0 = exp(log(10)*logy0);
	double chic1 = ChiCritical(x0, y0);


	char str[100];
	sprintf(str, "data/x_%d_y_%d", i, j);
	ofstream f1(str);
	f1 << chic1 << endl;
	f1.close();



//	if (0) {
//		for (int n = 0; n < loopy; n++) {
//			y0 = exp(log(10)*logy0);
//			logx0 = log10(initx);
//			for (int m = 0; m < loopx; m++) {
//				x0 = exp(log(10)*logx0);

//				chiCritical[m][n] = ChiCritical(x0, y0);
//				logx0 += hx;
//			}
//			logy0 += hy;

//		}
//	}
//	else {
//
//		double chiCritical1[5][5] = { {2.6154,2.5633,3.3945,3.8296,2.1542},
//		{2.5511,2.7729,4.3495,4.5542,3.4050},{2.6892,3.7866,4.5823,4.5523,4.2425},
//		{2.8240,4.1058,4.8172,4.529,4.8116},{3.571,4.1772,4.4678,4.5948,4.6152} };
//		memcpy(chiCritical, chiCritical1, sizeof(chiCritical1));

//	}





	


	return 0;




}




