
#define _CRT_SECURE_NO_WARNINGS


//#include<TString.h>
#include<TRandom3.h>
#include<TMinuit.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<algorithm>
using namespace std;

TRandom3 rnd(0);
int const bin = 5;
int const background = 100;
int const Nn = 1000;
double Ebins[] = { 10, 20, 30, 40, 50, 60 };
double L1 = 0.6, L2 = 1;
static int nExperiment[bin];
static int *nexp = &nExperiment[0];
double initx = 0.0001;
double xlower = -4.;
int nx = 200, ny = 200;
int const loopx = 10, loopy = 10;



double Mu(double x, double y, int binNumber)
{
	double E2 = Ebins[binNumber + 1];
	double E1 = Ebins[binNumber];


	auto Simpson = [&](int n) {//a is lower bound, b is upper bound

		double h = (E2 - E1) / n;
		double e = E1, sum0 = 0, sum = 0;
		for (int i = 0; i < n; i++)
		{
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
	int const N = 3000;
	double muNumber = x * Simpson(N) * 10000 / (L2 - L1) / (E2 - E1);// n = 1000, which is the number of integration bins

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


void DrawDeltaChiC(double(&ChiC)[loopx][loopy]) 
{
	TCanvas cvs;
	cvs.SetWindowSize(200, 300);

	int lbx = int(xlower);
	TH2D hist("", "", loopx, lbx, 0, loopy, 0, 3);
	for (int xi = 0; xi < loopx; ++xi) {
		for (int yi = 0; yi < loopy; ++yi) {
			//cout << xi << ";" << yi << endl;
			hist.SetBinContent(xi + 1, yi + 1, ChiC[xi][yi]);
		}
	}
	gStyle->SetOptStat(0);
	hist.Draw("TEXT");
	cvs.Print("ChiCTEXT2.png");
	hist.Draw("COLZ");
	cvs.Print("ChiCCOLZ2.png");
}


double Interpolate(double(&ChiC)[loopx][loopy], double x, double y) {
	//cout << x << "x" << y << "y" << endl;
	double xt = log10(x), yt = log10(y);
	double xid = (xt - xlower) / (-xlower) * (loopx - 1);
	double yid = yt / 3 * (loopy - 1);
	int xlow = floor(xid);
	int ylow = floor(yid);

	//cout << xlow << ";" << ylow << endl;
	if (xlow < 0) {
		xlow = 0;
		xid = xlow;
	}
	else  if (xlow >= loopx - 1) {
		xlow = loopx - 2;
		xid = xlow + 1;
	}
	if (ylow < 0) {
		ylow = 0;
		yid = ylow;
	}
	else  if (ylow >= loopy - 1) {
		ylow = loopy - 2;
		yid = ylow + 1;
	}

	int xup = xlow + 1;
	int yup = ylow + 1;
	double xlowWei = xup - xid;
	double ylowWei = yup - yid;
	//cout << xid << ";" << xup << ";" << xlow << endl;
	double chiCr = xlowWei*ylowWei*ChiC[xlow][ylow]
		+ xlowWei*(1 - ylowWei)*ChiC[xlow][yup]
		+ (1 - xlowWei)*ylowWei*ChiC[xup][ylow]
		+ (1 - xlowWei)*(1 - ylowWei)*ChiC[xup][yup];

	return chiCr;

}


void DrawInter(double(&ChiC)[loopx][loopy]) {

	TCanvas cvs;
	gStyle->SetOptStat(0);
	//cvs.SetWindowSize(200, 300);


	TH2D hist("", "", nx, xlower, 0, ny, 0, 3);
	//double xi1 = 0, yi1 = 0;
	for (int xi = 0; xi < nx; ++xi) {
		for (int yi = 0; yi < ny; ++yi) {
			//xi1 = xi;
			//yi1 = yi;
			double xcare = exp(log(10)*(1.*xi / nx * (-xlower) + xlower));
			double ycare = exp(log(10)*(1.*yi / ny * 3));
			//double chi = Chi(xcare, ycare);
			//cout << "xi" << xi << "yi" << yi << endl;

			double chic = Interpolate(ChiC, xcare, ycare);

			hist.SetBinContent(xi + 1, yi + 1, chic);
		}
	}
	//double levels[] = { -1, 0,10 };
	//hist.SetContour(3, levels);
	//hist.Draw("TEXT");
	//cvs.Print("ChiCInterpolateTEXT.png");
	hist.Draw("COLZ");
	cvs.Print("ChiCInterpolateCOLZ.png");

}


void DrawContour(double(&ChiC)[loopx][loopy]) {
	TCanvas cvs;
	gStyle->SetOptStat(0);
	//cvs.SetWindowSize(200, 300);


	TH2D hist("", "", nx, xlower, 0, ny, 0, 3);
	//cout << "try_0" << endl;
	for (int xi = 0; xi < nx; ++xi) {
		for (int yi = 0; yi < ny; ++yi) {
			//cout << xi << ";" << yi << endl;
			double xcare = exp(log(10)*(1.*xi / nx*(-xlower) + xlower));
			double ycare = exp(log(10)*(1.*yi / ny * 3));
			double chi = Chi(xcare, ycare);
			double chic = Interpolate(ChiC, xcare, ycare);
			//cout << chic << endl;
			double difference = chi - chic;
			//cout << "chic = " << chic << endl;
			//cout << "chi - chic = " << difference << endl;

			if (difference > 15) { difference = 15; }
			hist.SetBinContent(xi + 1, yi + 1, difference);
		}
	}
	double levels[] = { -10,0 };
	hist.SetContour(2, levels);
	//cout << "try" << endl;
	//hist.GetZaxis()->SetRangeUser(-4.61,0);
	hist.Draw("SAME CONTZ");
	//hist.Draw("colz");
	cvs.Print("ChiContourCOLZ.png");




}


int main()
{
	double chiCritical[loopx][loopy] = {0};
	ifstream f1;
	char str[100];
	double chic1;
	for (int m = 0; m < loopx; m++) {
		for (int n = 0; n < loopy; n++) {
			sprintf(str, "Chidata/x_%dy_%d", m, n);
			f1.open(str);
			f1 >> chic1;
			chiCritical[m][n] = chic1;
			f1.close();

		}

	}


	double xTrue = 0.006;
	double yTrue = 40;


	double muT[bin];

	for (int i = 0; i < bin; i++) {
		muT[i] = Mu(xTrue, yTrue, i);

		nExperiment[i] = Ngenerated(muT[i]);
		cout << nExperiment[i] << endl;
		//nExperiment[i] = 100;
	}

	DrawDeltaChiC(chiCritical);
	DrawInter(chiCritical);
	DrawContour(chiCritical);


	R__EXTERN TStyle  *gStyle;
	return 0;

}