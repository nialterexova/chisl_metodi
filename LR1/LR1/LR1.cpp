﻿#include <iostream>
#include <fstream>

using namespace std;

void analitic(double Tmax, double Xmax, double hx, double dt, int Nx, int Nt, double** u) {
	ofstream f;
	f.open("analytic.csv");	//файл для записи резуьтатов, посчитанных аналитически
	ofstream ErrorFile;
	ErrorFile.open("error.csv"); //файл для записи результатов погрешности решения

	double pi = 3.1416;
	double pi2 = pi * pi;
	double c[8];

	double** v = new double* [Nt + 1];
	for (int i = 0; i < Nt + 1; i++) {
		v[i] = new double[Nx];		//создание массива
	}

	for (int k = 1; k < 8; k++) {
		if (k % 2 == 0)  c[k] = 0;
		else c[k] = -100 / (pi2 * k * k);		//определение коэффициента Ск
	}

	for (int n = 0; n <= Nt; n++) {
		double sumErr = 0;
		for (int j = 0; j < Nx; j++) {
			v[n][j] = 25 / 2;
			for (int k = 1; k < 8; k++) v[n][j] += c[k] * exp(-pi2 * k * k * 1.14 * n * dt / (25*25)) * cos(pi * k * j * hx / 25);	//определение аналит решения
			f << v[n][j] << "; ";		//вывод в файл аналит решения точки 
			sumErr += (v[n][j] - u[n][j]) * (v[n][j] - u[n][j]);		//суммирование ошибки для данной итерации по времени
		}
		f << endl;
		ErrorFile << n * dt << "; " << sqrt(sumErr) << endl;		//вывод в файл ошибки для итерации по времени
	}
	f.close();
	ErrorFile.close();
}

int main()
{
	ofstream file;
	file.open("result.csv");

	double Xmax, Tmax, hx, dt;
	Xmax = 24; Tmax = 180; hx = 1; dt = 0.05;			//задание начальных условий

	int Nx = 1 + Xmax / hx;			//вычисление количества итераций
	int Nt = Tmax / dt;

	int j, n = 0;
	double** u = new double* [Nt + 1];
	for (int i = 0; i < Nt + 1; i++) {
		u[i] = new double[Nx];		//выделение памяти под массив 
	}

	for (j = 0; j < Nx; j++) {
		u[0][j] = hx * j;		//начальные условия
	}

	while (n < Nt)
	{
		for (j = 1; j < Nx - 1; j++) {
			u[n + 1][j] = u[n][j] + (1.14 * dt * (u[n][j + 1] - 2 * u[n][j] + u[n][j - 1])) / (hx * hx);	//рекуррентное соотношение
		}
		u[n + 1][0] = u[n + 1][1];					//определение первой и последней точки по Х из граничных условий
		u[n + 1][Nx - 1] = u[n + 1][Nx - 2];
		n++;
	}

	for (n = 0; n <= Nt; n++) {
		for (j = 0; j < Nx; j++) {
			file << u[n][j] << "; ";	//вывод результата 
		}
		file << endl;
	}

	analitic(Tmax, Xmax, hx, dt, Nx, Nt, u);	//Запуск аналитического решения

	for (int i = 0; i <= Nt; i++) {
		delete[]u[i];
	}

	file.close();
}
