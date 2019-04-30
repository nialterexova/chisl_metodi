
#include <iostream>
#include <fstream>

using namespace std;

void analit(double Tmax, double Xmax, double hx, double ht, int Nx, int Nt, double** U) {
	ofstream f;
	f.open("analyt.csv");
	ofstream Err;
	Err.open("error.csv");

	double pi = 3.1416;
	double pi2 = pi * pi;
	double c[8];
	double** v = new double* [Nt + 1];
	for (int i = 0; i < Nt + 1; i++) {
		v[i] = new double[Nx];
	}
	for (int k = 1; k < 8; k++) {
		if (k % 2 == 0)  c[k] = 0;
		else c[k] = 60 * 2 / (pi2 * k * k);
	}
	for (int n = 0; n <= Nt; n++) {
		double sumErr = 0;
		for (int j = 0; j < Nx; j++) {
			v[n][j] = 30/2;
			for (int k=1; k < 8; k++) v[n][j] += c[k] * exp(-pi2 * k * k * 10 * n * ht / 900) * cos(pi * k * j * hx / 30);
			f << v[n][j] << "; ";
			sumErr += (v[n][j] - U[n][j]) * (v[n][j] - U[n][j]);
		}
		f << endl;
		Err << n * ht << "; " << sqrt(sumErr) << endl;
	}
	f.close();
	Err.close();
}

int main()
{
	ofstream fout;
	fout.open("lr1.csv");

	double Tmax, Xmax = 30;
	double hx = 0, ht = 0;
	/*cout << "Enter Tmax" << endl;
	cin >> Tmax;
    cout << "Enter hx and ht" << endl;
	cin >> hx >> ht;*/
	Tmax = 150; hx = 1; ht = 0.05;

	
	int Nx = 1 + Xmax / hx;
	int Nt = Tmax / ht;
	fout << "Nx = " << Nx << endl;
	fout << "Nt = " << Nt << endl;

	int j, n = 0;
	double** U = new double* [Nt+1];
	for (int i = 0; i < Nt+1; i++) {
		U[i] = new double[Nx];
	}

	for (j = 0; j < Nx; j++) {
		U[0][j] = 30 - hx * j ;
	}

	while (n<Nt)
	{
		for (j = 1; j < Nx-1; j++) {
			U[n + 1][j] = U[n][j] + (10 * ht * (U[n][j + 1] - 2 * U[n][j] + U[n][j - 1])) / (hx * hx);
		}
		U[n + 1][0] = U[n + 1][1];
		U[n + 1][Nx-1] = U[n + 1][Nx-2];
		n++;
	}

	for (n = 0; n <= Nt; n++) {
		if (n == 0) {
			fout << "n/j; ";
			for (j = 0; j < Nx; j++) fout << hx * j << "; ";
			fout << endl;
		}
		fout << n << "(" << n*ht <<"); ";
		for (j = 0; j < Nx; j++) {
			fout << U[n][j] << "; ";
		}
		fout << endl;
	}

	//--------analit-------------

	analit(Tmax, Xmax, hx, ht, Nx, Nt, U);

	//--------analit-------------

	for (int i = 0; i <= Nt; i++) {
		delete[]U[i];
	}

	fout.close();
}
