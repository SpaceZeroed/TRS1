// TRS1.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
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
	inline double boundcondf(double x)
	{
		return exp(x) * (x * x + x - 3);
	}
	inline double p(double x)
	{
		return -2.;
	}
	inline double q(double x)
	{
		return 0.;
	}
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

vector <double> Diag3Prog(vector<double> a, vector<double> b, vector<double> c , vector<double> f , int N)
{
	vector <double> y(N, { 0 });
	vector <double> xi(N + 1, { 0. });
	vector <double> eta(N + 1, { 0. });
	// a - b + c = d
	for (int i = 0; i < N; i++) // прямой ход 
	{
		xi[i + 1] = -c[i] / (a[i] * xi[i] + b[i]);
		eta[i + 1] = (f[i] - eta[i] * a[i]) / (a[i] * xi[i] + b[i]);
	}
	y[N - 1] = eta[N];
	for (int i = N - 1; i > 0; i--) // обратный ход 
	{
		y[i - 1] = xi[i] * y[i] + eta[i];
	}
	return y;
}
vector<pair<double, double>> FiniteDifferenceMethod(double N, double start_point, double end_point)
{
	double alpha_n = 0., beta_n = 1., gamma = 1., delta = 1., A = 2., B = 2.;
	double h = (end_point - start_point) / N;

	vector<pair<double, double>> tx;
	vector <double> a(N, { 0 });
	vector <double> b(N, { 0 });
	vector <double> c(N, { 0 });
	vector <double> d(N, { 0 });
	
	a[0] = 0.;
	b[0] = alpha_n - beta_n / h;
	c[0] = beta_n / h;
	d[0] = A;
	for (int i = 1; i < N - 1; i++)
	{
		a[i] = 1. / (h * h) - p(i * h) / (2. * h);
		b[i] = q(i * h) - 2. / (h * h);
		c[i] = 1. / (h * h) + p(i * h) / (2 * h);
		d[i] = boundcondf(i * h);
	}
	a[N - 1] = -delta / h;
	b[N - 1] = gamma + delta / h;
	c[N - 1] = 0.;
	d[N - 1] = B;

	vector<double> y = Diag3Prog(a, b, c, d, N);
	for (int i = 0; i < N; i++)
		tx.push_back(make_pair(start_point + h * i, y[i] ));
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
	ofstream Ex5("Ex5_xuv.txt");
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
	double h = 0.001;

	//tvx = EilerMeth(h);
	//Outfile("EilerMethod", tvx, h);

	//tvx = AdamsMethod(h);
	//Outfile("AdamsMethod", tvx, h);

	//tvx = RungeKuttMethod(h);
	//Outfile("RKMethod", tvx, h);
	//Outcmd(tvx);

	tx = FiniteDifferenceMethod(250, 0, 1);
	OutfileTX("RKMethod", tx, h);

	//Ex5();
	return 0;
}
