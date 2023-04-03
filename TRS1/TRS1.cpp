// TRS1.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
// 
#include <iostream> 
#include <math.h> 
#include <vector> 
#include <tuple> 
#include <fstream> 
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
	double bcf(double x)
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
		double xn, vn, tn, tn_1,vn_1, xn_1;

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
	while (t<10)
	{
		double tempV, tempX;
		t += h;
		tempX = Xn + alpha*(k1 + 2 * k2 + 2 * k3 + k4) / 6;
		tempV = Vn + betta*(m1 + 2 * m2 + 2 * m3 + m4) / 6;
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

vector<tuple<double, double, double>> FiniteDifferenceMethod(double h,double a,double b)
{
	int n = int((b - a) / h);

}
int main()
{

	vector<tuple<double, double, double>> tvx;
	double h = 0.001;

	//tvx = EilerMeth(h);
	//Outfile("EilerMethod", tvx, h);

	//tvx = AdamsMethod(h);
	//Outfile("AdamsMethod", tvx, h);

	tvx = RungeKuttMethod(h);
	Outfile("RKMethod", tvx, h);
	//Outcmd(tvx);
	return 0;
}
