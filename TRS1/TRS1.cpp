﻿// TRS1.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
// 
#include <iostream> 
#include <math.h> 
#include <vector> 
#include <tuple> 
#include <fstream>
#include <iomanip>
#include <string> 
namespace var9
{
	// Cauchy's problem data
	double m = 6.65 * pow(10, -27);
	double x0 = 0;
	double v0 = 1;
	double U(double x)
	{
		return(6 * pow(10, -17) * x * x * sinh(x * x));
	}
	double f(double x)
	{
		return(x * sinh(x * x) + 2 * pow(x, 3) * cosh(x * x));
	}
	double alpha = 1;
	double betta = -12. / 665;
	// Boundary value's problem data
	double a0 = 0;
	double a1 = 1;
	double b0 = 1;
	double b1 = 1;
	double p = -2;
	double q = 0;
	double boundcondf(double x)
	{
		return exp(x) * (x * x + x - 3);
	}
	double A = 2;
	double B = 2;
}
using namespace std;
using namespace var9;

vector<tuple<double, double, double>> EilerMeth(double h)
{

	double t = 0;
	vector<tuple<double, double, double>> tvx;
	tvx.push_back(make_tuple(t, v0, x0));
	while (t < 10)
	{
		double tempV, tempX;
		double Xn, Vn, Tn;
		tie(Tn, Vn, Xn) = tvx[tvx.size() - 1];
		tempV = Vn + h * betta * f(Xn);
		tempX = Xn + h * alpha * Vn;
		t += h;
		tvx.push_back(make_tuple(t, tempV, tempX));
	}
	return tvx;
}
void Outcmd(vector<tuple<double, double, double>> tvx)
{
	cout << "table of contents: t v x columns\n";
	for (int i = 0; i < tvx.size(); i++)
	{
		double tempT, tempX, tempV;
		tie(tempT, tempV, tempX) = tvx[i];
		cout << tempT << " " << tempV << " " << tempX << endl;
	}
}
void Outfile(string Name, vector<tuple<double, double, double>> tvx, double h)
{
	string filename = "Ex1His_" + Name + to_string(float(h)) + ".txt";
	string filename2 = "Tvex1His_" + Name + to_string(float(h)) + ".txt";
	string filename3 = "Txex1His_" + Name + to_string(float(h)) + ".txt";
	ofstream output(filename);
	ofstream output2(filename2);
	ofstream output3(filename3);
	for (int i = 0; i < tvx.size(); i++)
	{
		double tempT, tempX, tempV;
		tie(tempT, tempV, tempX) = tvx[i];
		output << tempT << " " << tempV << " " << tempX << endl;
		output2 << tempT << " " << tempV << endl;
		output3 << tempT << " " << tempX << endl;

	}
	output.close();
	output2.close();
	output3.close();
}
void OutfileTX(string Name, vector<pair<double, double>> tx, double h)
{
	string filename = "Ex4H_" + Name + to_string(float(h)) + ".txt";
	ofstream output(filename);
	for (int i = 0; i < tx.size(); i++)
	{
		output << tx[i].first << " " << tx[i].second << endl;
	}
	output.close();
}
vector<tuple<double, double, double>> EilerRungeMeth(double h)
{
	double q1 = 0.5;
	double h1 = 0.001;
	double h2 = h1 * q1;
	double p1 = 1;
	vector<tuple<double, double, double>> tvx1 = EilerMeth(h1);
	vector<tuple<double, double, double>> tvx2 = EilerMeth(h2);
	vector<tuple<double, double, double>> tvx;
	for (int i = 1; i < tvx1.size() - 1; i++)
	{
		double X1, X2, V1, V2, T2, T1;
		tie(T1, V1, X1) = tvx1[i];
		tie(T2, V2, X2) = tvx2[i * 2];
		tvx.push_back(make_tuple(T1, V1 + (V1 - V2) / (2 - 1), 2 * X1 - X2));
	}
	return tvx;
}
vector<tuple<double, double, double>> AdamsMethod(double h)
{

	vector<tuple<double, double, double>> tvx;
	double t = 0;
	tvx.push_back(make_tuple(t, v0, x0));
	tvx.push_back(make_tuple(t + h, v0 + h * betta * f(x0), x0 + h * alpha * v0));
	while (t < 10)
	{
		double tempV, tempX;
		double xn, vn, tn, tn_1, vn_1, xn_1;

		tie(tn, vn, xn) = tvx[tvx.size() - 1];
		tie(tn_1, vn_1, xn_1) = tvx[tvx.size() - 2];

		tempV = vn + h * betta * (3 * f(xn) - f(xn_1)) / 2;
		tempX = xn + alpha * (3 * vn - vn_1) * h / 2;
		t += h;
		tvx.push_back(make_tuple(t, tempV, tempX));
	}
	return tvx;
}
vector<tuple<double, double, double>> RungeKutty(double h)
{
	using namespace var9;
	vector<tuple<double, double, double>> tvx;
	double t = 0;
	tvx.push_back(make_tuple(t, v0, x0));
	while (t < 10)
	{
		double tempV, tempX;
		double xn, vn, tn, k_11, k_12, k_13, k_14, k_21, k_22, k_23, k_24;
		tie(tn, vn, xn) = tvx[tvx.size() - 1];
		k_11 = vn; k_12 = vn + k_11 * h; k_13 = vn + k_12 * h; k_14 = vn + k_13 * h;
		k_21 = alpha * f(xn); k_22 = alpha * f(xn + k_21 * h / 2); k_23 = alpha * f(xn + k_22 * h / 2); k_24 = alpha * f(xn + k_23 * h);
		tempX = xn + (k_11 + 2 * k_12 + 2 * k_13 + k_14) * h / 6;
		tempV = vn + (k_21 + 2 * k_22 + 2 * k_23 + k_24) * h / 6;
		t += h;
		tvx.push_back(make_tuple(t, tempV, tempX));
	}
	return tvx;
}
vector<tuple<double, double, double>> RungeKuttMethod(double h)
{
	vector<tuple<double, double, double>> tvx;
	double t = 0;
	tvx.push_back(make_tuple(t, v0, x0));
	double Xn, Vn, Tn;
	tie(Tn, Vn, Xn) = tvx[tvx.size() - 1];
	// kn - coeffs for X where Vn is derrative, mn - coeffs for V where f(xn) is derrative 
	double k1 = Vn * h;
	double m1 = f(Xn) * h;
	double k2 = (Vn + m1 / 2) * h;
	double m2 = f(Xn + k1 / 2) * h;
	double k3 = (Vn + m2 / 2) * h;
	double m3 = f(Xn + k2 / 2) * h;
	double k4 = (Vn + m3) * h;
	double m4 = f(Xn + k3) * h;
	while (t < 10)
	{
		double tempV, tempX;
		t += h;
		tempX = Xn + alpha * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
		tempV = Vn + betta * (m1 + 2 * m2 + 2 * m3 + m4) / 6;
		tvx.push_back(make_tuple(t, tempV, tempX));
		tie(Tn, Vn, Xn) = tvx[tvx.size() - 1];
		k1 = Vn * h;
		m1 = f(Xn) * h;
		k2 = (Vn + m1 / 2) * h;
		m2 = f(Xn + k1 / 2) * h;
		k3 = (Vn + m2 / 2) * h;
		m3 = f(Xn + k2 / 2) * h;
		k4 = (Vn + m3) * h;
		m4 = f(Xn + k3) * h;
	}
	return tvx;
}
// functions for boundary problem
void PrintMatrix(vector<vector<double>> Matrix)
{
	cout << fixed << std::setprecision(4);
	cout << "-------------------------------------------------------------" << endl;
	for (int i = 0; i < Matrix.size(); i++)
	{
		for (int j = 0; j < Matrix[i].size(); j++)
		{
			cout << setw(5) << Matrix[i][j] << "  ";
		}
		cout << endl;
	}
	cout << "-------------------------------------------------------------" << endl;
}

vector <double> Diag3Prog(vector<vector<double>> matrix, vector<double> f)
{
	int n = f.size();

	vector <double> x(n + 1, { 0 });
	vector <double> psi(n + 1, { 0 });
	vector <double> ksi(n + 1, { 0 });

	double b = -matrix[0][0]; double c = matrix[0][1];
	psi[1] = -c / (0 * psi[0] - b);
	ksi[1] = (f[0] - 0 * ksi[0]) / (0 * psi[0] - b);

	for (int i = 1; i < n - 2; i++) // прямой ход 
	{
		double a = matrix[i][i - 1]; double b = -matrix[i][i]; double c = matrix[i][i + 1];
		psi[i + 1] = -c / (a * psi[i] - b);
		ksi[i + 1] = (f[i] - a * ksi[i]) / (a * psi[i] - b);
	}
	for (int i = n; i > 0; i--) // обратный ход 
	{
		x[i - 1] = psi[i] * x[i] + ksi[i];
	}
	return x;
}
vector<pair<double, double>> FiniteDifferenceMethod(double h, double a, double b)
{
	int n_big = int((b - a) / h);
	vector<pair<double, double>> tx;
	vector<vector<double>> matrix_prog;
	matrix_prog.resize(n_big);
	for (int i = 0; i < n_big; i++)
		matrix_prog[i].resize(n_big);

	matrix_prog[0][0] = -1 / h;
	matrix_prog[0][1] = 1 / h;
	for (int i = 1; i <= n_big - 2; i++)
	{
		matrix_prog[i][i - 1] = 2 - h * p;
		matrix_prog[i][i] = -4 + q * 2 * h * h;
		matrix_prog[i][i + 1] = 2 + h * p;
	}
	matrix_prog[n_big - 1][n_big - 1] = 1 / h + 1;
	matrix_prog[n_big - 1][n_big - 2] = -1 / h;

	vector<double> f;
	f.resize(n_big);
	f[0] = 2; f[n_big - 1] = 2;
	for (int i = 1; i < n_big - 1; i++)
	{
		f[i] = boundcondf(a + i * h);
	}
	vector<double> u = Diag3Prog(matrix_prog, f);
	for (int i = 0; i < n_big; i++)
		tx.push_back(make_pair(a + h * i, u[i]));
	//PrintMatrix(matrix_prog);
	return tx;
}

double F(double v)
{
	return v;
}

double G( double x,double v)
{
	return pow(exp(x*x + x - 3), x) + v;
}

double RK4(double p, double n)
{
	using namespace var9;
	double k0, k1, k2, k3, m0, m1, m2, m3;
	double h = 1. / n;
	double x_0, v_0, u_n, v_n, u_0;
	
	x_0 = 0;
	u_0 = p;
	v_0 = 2;
	k0 = F(v_0) * h;
	m0 = G(x_0, v_0) * h; 
	k1 = F(v_0 + m0/2) * h;
	m1 = G(x_0 + h / 2, v_0 + m0 / 2) * h;
	k2 = F(v_0 + m1 / 2) * h;
	m2 = G(x_0 + h / 2, v_0 + m1 / 2) * h;
	k3 = F(v_0 + m2) * h;
	m3 = G(x_0 + h, v_0 + m2) * h;
	for (int i = 0; i < n; i++)
	{
		x_0 = i * h;
		u_n = u_0 + (k0 + k1 + k2 + k3) / 6;
		v_n = v_0 + (m0 + m1 + m2 + m3) / 6;
		u_0 = u_n;
		v_0 = v_n;
		k0 = F(v_0) * h;
		m0 = G(x_0, v_0) * h;
		k1 = F(v_0 + m0 / 2) * h;
		m1 = G(x_0 + h / 2, v_0 + m0 / 2) * h;
		k2 = F(v_0 + m1 / 2) * h;
		m2 = G(x_0 + h / 2, v_0 + m1 / 2) * h;
		k3 = F(v_0 + m2) * h;
		m3 = G(x_0 + h, v_0 + m2) * h;
	}
	
	return (u_0 + v_0);
}
double find_p(double n)
{
	double p1 = 10;
	double p2 = -10.;
	double p = (p2 - p1) / 2.;
	while (RK4(p1, n) - 2 > 1e-10)
	{
		cout << p;
		p = (p2 + p1) / 2.;
		if ((RK4(p1, n) - 2) * (RK4(p, n) - 2) < 0)
		{
			p2 = p;
		}
		else if ((RK4(p, n) - 2 ) * (RK4(p2, n) - 2 ) < 0)
		{
			p1 = p;
		}
	}
	return p1;
}
void Ex5()
{
	ofstream Ex5("Ex5_xuv");
	vector<pair<double, double>> tx;
	double p = find_p(1000);
	double  n = 1000;
	using namespace var9;
	double k0, k1, k2, k3, m0, m1, m2, m3;
	double h = 1. / n;
	double x_0, v_0, u_n, v_n, u_0;
	x_0 = 0;
	u_0 = p;
	v_0 = 2;
	Ex5 << x_0 << " " << u_0 << endl;
	k0 = F(v_0) * h;
	m0 = G(x_0, v_0) * h;
	k1 = F(v_0 + m0 / 2) * h;
	m1 = G(x_0 + h / 2, v_0 + m0 / 2) * h;
	k2 = F(v_0 + m1 / 2) * h;
	m2 = G(x_0 + h / 2, v_0 + m1 / 2) * h;
	k3 = F(v_0 + m2) * h;
	m3 = G(x_0 + h, v_0 + m2) * h;
	for (int i = 0; i < n; i++)
	{
		x_0 = i * h;
		u_n = u_0 + (k0 + k1 + k2 + k3) / 6;
		v_n = v_0 + (m0 + m1 + m2 + m3) / 6;
		u_0 = u_n;
		v_0 = v_n;
		Ex5 << x_0 << " " << u_0 << endl;
		k0 = F(v_0) * h;
		m0 = G(x_0, v_0) * h;
		k1 = F(v_0 + m0 / 2) * h;
		m1 = G(x_0 + h / 2, v_0 + m0 / 2) * h;
		k2 = F(v_0 + m1 / 2) * h;
		m2 = G(x_0 + h / 2, v_0 + m1 / 2) * h;
		k3 = F(v_0 + m2) * h;
		m3 = G(x_0 + h, v_0 + m2) * h;
	}
	Ex5.close();
}
int main()
{

	vector<tuple<double, double, double>> tvx;
	vector<pair<double, double>> tx;
	double h = 0.01;

	//tvx = EilerMeth(h);
	//Outfile("EilerMethod", tvx, h);

	//tvx = AdamsMethod(h);
	//Outfile("AdamsMethod", tvx, h);

	//tvx = RungeKuttMethod(h);
	//Outfile("RKMethod", tvx, h);
	//Outcmd(tvx);

	//tx = FiniteDifferenceMethod(h, 0, 1);
	//OutfileTX("RKMethod", tx, h);

	Ex5();
	return 0;
}
