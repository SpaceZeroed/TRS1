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
}
using namespace std;
//double dxdt(double v)
//{
//    return v;
//}


vector<tuple<double, double, double>> EilerMeth(double h)
{
	using namespace var9;
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
void Outfile(vector<tuple<double, double, double>> tvx, double h)
{
	string filename = "Ex1His_" + to_string(float(h)) + ".txt";
	string filename2 = "Tvex1His_" + to_string(float(h)) + ".txt";
	string filename3 = "Txex1His_" + to_string(float(h)) + ".txt";
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
	double q = 0.5;
	double h1 = 0.001;
	double h2 = h1 * q;
	double p = 1;
	vector<tuple<double, double, double>> tvx1 = EilerMeth(h1);
	vector<tuple<double, double, double>> tvx2 = EilerMeth(h2);
	vector<tuple<double, double, double>> tvx;
	for (int i = 1; i < tvx1.size() - 1; i++)
	{
		double X1, X2, V1, V2, T2, T1;
		tie(T1, V1, X1) = tvx1[i];
		tie(T2, V2, X2) = tvx2[i * 2];
		tvx.push_back(make_tuple(T1,V1 + (V1 - V2) / (2 - 1), 2 * X1 - X2));
	}

}
vector<tuple<double, double, double>> AdamsMeth(double h)
{
	double h = 0.001;
	vector<tuple<double, double, double>> tvx = EilerMeth(h);

}

int main()
{
	using namespace var9;
	vector<tuple<double, double, double>> tvx;
	double h = 0.001;
	tvx = EilerMeth(h);
	//Outcmd(tvx);
	Outfile(tvx, h);
	return 0;
}
// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
