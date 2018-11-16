// HMM.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"

//#define DEBUG_MODE

#ifdef DEBUG_MODE
#define DEBUG(x) cout<<#x<<"\t: "<<x<<std::endl
#define DEBUGV(x) cout<<#x<<"\t: {"; for(auto _x : x) cout<<'\t'<<_x; cout<<"}"<<std::endl
#define DEBUGA(x,s,e) cout<<#x<<'['<<s<<':'<<e<<"]\t: {";for(int _i=s;_i<e;)cout<<x[_i++]<<'\t';cout<<"\b }"<<std::endl
#else
#define DEBUG(x)
#define DEBUGV(x)
#define DEBUGA(x,s,e)
#endif // DEBUG_MODE

using vec_int = std::vector<int>;
using vec_dbl = std::vector<double>;
using matrix_dbl = std::vector<vec_dbl>;
using matrix_int = std::vector<vec_int>;
using std::string;
using std::cout;
using std::cin;

void read_mat(matrix_dbl &m, string infile, int breadth)
{
	std::ifstream fin(infile);
	int inp;
	if (!fin)
		throw std::runtime_error(infile + " not found!");
	while (true)
	{
		fin >> inp;
		if (fin.eof())
			return;
		m.push_back(vec_dbl(breadth));
		for (int i = 0; i < breadth; i++)
			fin >> m.back()[i];
	}
}

void read_obs(vec_int &m, string infile)
{
	std::ifstream fin(infile);
	int inp;
	if (!fin)
		throw std::runtime_error(infile + " not found!");
	while (m.size() < 85)
	{
		fin >> inp;
		if (fin.eof())
			return;
		m.push_back(inp);
	}
}

void read_pi(vec_dbl &m, string infile)
{
	std::ifstream fin(infile);
	double inp;
	if (!fin)
		throw std::runtime_error(infile + " not found!");
	while (true)
	{
		fin >> inp;
		if (fin.eof())
			return;
		m.push_back(inp);
	}
}

int add(int a, int b)
{
	return a + b;
}

int main()
{
	//int line = 72; DEBUG(line);
	vec_dbl pi;
	matrix_dbl a, b, alpha, beta, gamma, delta;
	matrix_int psi;
	vec_int o, q;
	try
	{
		read_mat(a, (string)inpath + "a", NUM_STATES);
		read_mat(b, (string)inpath + "b", M);
		read_pi(pi, (string)inpath + "pi");
		read_obs(o, (string)inpath + "obs");
	}
	catch (std::runtime_error e)
	{
		cout << e.what() << " Aborting...\n";
		return -1;
	}

	for (auto ia : a)
	{
		DEBUGV(ia);
	}
	for (auto ib : b)
	{
		DEBUGV(ib);
	}
	DEBUGV(pi);
		
	alpha.resize(o.size(), vec_dbl(NUM_STATES, 0));
	beta.resize(o.size(), vec_dbl(NUM_STATES, 0));
	gamma.resize(o.size(), vec_dbl(NUM_STATES, 0));
	delta.resize(o.size(), vec_dbl(NUM_STATES, 0));
	psi.resize(o.size(), vec_int(NUM_STATES, 0));
	q.resize(o.size());
	for (int i = 0; i < NUM_STATES; i++)
	{
		alpha[0][i] = pi[i] * b[i][o[0]];
	}
	for (int t = 1; t < o.size(); t++)
	{
		for (int i = 0; i < NUM_STATES; i++)
		{
			for (int j = 0; j < NUM_STATES; j++)
			{
				alpha[t][i] += alpha[t - 1][j] * a[j][i];
			}
			alpha[t][i] *= b[i][o[t]];
		}
	}
	for (int i = 0; i < NUM_STATES; i++)
		beta.back()[i] = 1;
	for (int t = o.size() - 2; t >= 0; t--)
	{
		for (int i = 0; i < NUM_STATES; i++)
		{
			for (int j = 0; j < NUM_STATES; j++)
			{
				beta[t][i] += beta[t + 1][j] * a[i][j] * b[j][o[t + 1]];
			}
		}
	}
	double sum;
	for (int t = 0; t < o.size(); t++)
	{
		sum = 0;
		for (int i = 0; i < NUM_STATES; i++)
		{
			sum += alpha[t][i] * beta[t][i];
		}
		for (int i = 0; i < NUM_STATES; i++)
		{
			gamma[t][i] = alpha[t][i] * beta[t][i] / sum;
		}
	}

	int maxp;
	double maxd, tempd;
	for (int i = 0; i < NUM_STATES; i++)
	{
		delta[0][i] = pi[i] * b[i][o[0]];
		psi[0][i] = 0;
	}
	for (int t = 1; t < o.size(); t++)
	{
		for (int j = 0; j < NUM_STATES; j++)
		{
			maxd = -1;
			maxp = 0;
			for (int i = 0; i < NUM_STATES; i++)
			{
				tempd = delta[t - 1][i] * a[i][j];
				if (maxd < tempd)
				{
					maxp = i;
					maxd = tempd;
				}
			}
			delta[t][j] = maxd * b[j][o[t]];
			psi[t][j] = maxp;
		}
	}
	maxp = 0;
	for (int i = 1; i < NUM_STATES; i++)
		if (delta.back()[i] > delta.back()[maxp]) maxp = i;
	double pstar = delta.back()[maxp];
	q.back() = maxp;
	for (int t = o.size() - 2; t >= 0; t--)
		q[t] = psi[t + 1][q[t + 1]];
	cout <<"pstar: "<< pstar << std::endl<< "qstar: ";
	for (auto iq : q)
		cout << iq << ' ';
	cout << std::endl;
	return 0;
}
