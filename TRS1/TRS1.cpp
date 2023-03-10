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
    double m = 6.65*pow(10,-27);
    double x0 = 0;
    double v0 = 1;
    double U(double x)
    {
        return(6*pow(10,-17) *x*x* sinh(x * x));
    }
    double dudx(double x)
    {
        return(6*pow(10,-17) * (2 * x*sinh(x*x)+2*pow(x,3)*cosh(x*x)));
    }
    
}
using namespace std;
double dxdt(double v)
{
    return v;
}

double dvdt(double x) // dont understand why i need x
{
    using namespace var9;
    return(-dudx(x) / m);
}
vector<tuple<double, double, double>> EilerFor2(double h)
{
    using namespace var9;
    double t = 0;
    // t is first, x is second 
    vector<pair<double, double>> tx;
    // t is first, v is second
    vector<pair<double, double>> tv;
    // finding v(t) first, then x(t)
    tv.push_back(make_pair(t, v0));
    tx.push_back(make_pair(t, x0));
    t += h*pow(10,-6);
    while (t < 1)
    {
        double tempV = 0;
        double tempX = 0;
        tempV = tv[tv.size() - 1].second + h * dvdt(tx[tx.size() - 1].second);
        tempX = tx[tx.size() - 1].second + h * dxdt(tv[tv.size() - 1].second);
        t += h;
        tv.push_back(make_pair(t, tempV));
        tx.push_back(make_pair(t, tempX));
    }
    // tuple for tvx
    vector<tuple<double, double, double>> tvx;
    for (int it = 0; it < tv.size(); it++)
    {
        tvx.push_back(make_tuple(tv[it].first, tv[it].second, tx[it].second));
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
        //tempV=tempV*pow(10,-40);
        //tempX=tempX*pow(10,-40);
        cout << tempT << " " << tempV << " " << tempX << endl;
    }
}
void Outfile(vector<tuple<double, double, double>> tvx, double h)
{
    h = h * 100;
    string filename = "Ex1His"+to_string(int(h))+".txt";
    ofstream output(filename);
    for (int i = 0; i < tvx.size(); i++)
    {
        double tempT, tempX, tempV;
        tie(tempT, tempV, tempX) = tvx[i];
        output << tempT << " " << tempV << " " << tempX << endl;

    }
    output.close();
}
int main()
{
    using namespace var9;
    vector<tuple<double, double, double>> tvx;
    double h = 0.1;
    tvx = EilerFor2(h);
    Outcmd(tvx);
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
