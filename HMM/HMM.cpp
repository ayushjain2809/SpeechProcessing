// HMM.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"

//#define DEBUG_MODE

#ifdef DEBUG_MODE
#define DEBUG(x) cout<<#x<<"\t: "<<x<<std::endl
#define DEBUGV(x) cout<<#x<<"\t: {"; for(auto _x : x) cout<<'\t'<<_x; cout<<"}"<<std::endl
#define DEBUGM(x) cout<<#x<<"\t: {"; for(auto _xr : x)  { cout<<std::endl; for(auto _x : _xr) cout<<'\t'<<_x; } cout<<"}"<<std::endl
#define DEBUGA(x,s,e) cout<<#x<<'['<<s<<':'<<e<<"]\t: {";for(int _i=s;_i<e;)cout<<x[_i++]<<'\t';cout<<"\b }"<<std::endl
#else
#define DEBUG(x)
#define DEBUGV(x)
#define DEBUGM(x)
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
	while (true)
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

void write_output(vec_int &v, string outfile)
{
	std::ofstream fout(outfile);
	for (auto x : v)
	{
		fout << x << ' ';
	}
	fout.close();
}

void write_output(vec_dbl &v, string outfile)
{
	std::ofstream fout(outfile);
	for (auto x : v)
	{
		fout << x << ' ';
	}
	fout.close();
}

void write_output(matrix_dbl &v, string outfile)
{
	std::ofstream fout(outfile);
	for (auto xr : v)
	{
		for (auto x : xr)
			fout << x << ' ';
		cout << '\n';
	}
	fout.close();
}

int prepare()
{
	//int line = 72; DEBUG(line);

	makeUniverse();

	matrix_dbl C = readUniverse();
	char dig[] = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' };
	vec_int LBGresult = LBGClassifier(C);
	write_output(LBGresult, (string)outpath + "observations");

	DEBUGV(LBGresult);
	std::ofstream fout(OutputFile);

	return 0;
}

double prob(vec_int &o, matrix_dbl &a, matrix_dbl &b, vec_dbl &pi)
{
	int T = o.size();
	matrix_dbl alpha, beta, gamma, delta;//, a_init, b_init, a_avg, b_avg, a_avg_prev, b_avg_prev;
	matrix_int psi;
	alpha.resize(T, vec_dbl(NUM_STATES, 0));
	beta.resize(T, vec_dbl(NUM_STATES, 0));
	gamma.resize(T, vec_dbl(NUM_STATES, 0));
	delta.resize(T, vec_dbl(NUM_STATES, 0));
	psi.resize(T, vec_int(NUM_STATES, 0));
	for (int i = 0; i < NUM_STATES; i++)
	{
		alpha[0][i] = pi[i] * b[i][o[0]];
	}
	for (int t = 1; t < T; t++)
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
	for (int t = T - 2; t >= 0; t--)
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

	int maxp;
	double maxd, tempd;
	for (int i = 0; i < NUM_STATES; i++)
	{
		delta[0][i] = pi[i] * b[i][o[0]];
		psi[0][i] = 0;
	}
	for (int t = 1; t < T; t++)
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
	return /*pstar =*/ delta.back()[maxp];
}

int train()
{
	vec_dbl pi, pi_init, pi_avg, pi_avg_prev;
	matrix_dbl a, b, alpha, beta, gamma, delta, a_init, b_init, a_avg, b_avg, a_avg_prev, b_avg_prev;
	matrix_int psi;
	vec_int o, q, LBGresult;
	try
	{
		read_mat(a_init, (string)inpath + "a", NUM_STATES);
		read_mat(b_init, (string)inpath + "b", M);
		read_pi(pi_init, (string)inpath + "pi");
		read_obs(LBGresult, (string)inpath + "observations");
	}
	catch (std::runtime_error e)
	{
		cout << e.what() << " Aborting...\n";
		return -1;
	}
	double pstar_prev;

	for (int d = 0; d < 10; d++)
	{
		a_avg = a_init;
		b_avg = b_init;
		pi_avg = pi_init;
		for (int train_iter = 0; train_iter < 3; train_iter++)
		{
			a_avg_prev = a_avg;
			b_avg_prev = b_avg;
			pi_avg_prev = pi_avg;
			for (int i = 0; i < NUM_STATES; i++)
			{
				for (int j = 0; j < NUM_STATES; j++)
				{
					a_avg[i][j] = 0;
				}
			}
			for (int i = 0; i < NUM_STATES; i++)
			{
				for (int j = 0; j < M; j++)
				{
					b_avg[i][j] = 0;
				}
			}
			for (int i = 0; i < NUM_STATES; i++)
			{
				pi_avg[i] = 0;
			}
			for (int file_num = 0; file_num < 6; file_num++)
			{
				cout << d << ' ' << file_num + 1 << '\n';
				a = a_avg_prev;
				b = b_avg_prev;
				pi = pi_avg_prev;
				o = vec_int(LBGresult.begin() + (d*num_of_rec + file_num)*no_of_itern, LBGresult.begin() + (d*num_of_rec + file_num + 1) * no_of_itern);
				int T = o.size();
				pstar_prev = 0;
				for (int iter = 0; iter < 25; iter++)
				{
					/*DEBUGM(a)
					DEBUGM(b)
					DEBUGV(pi);*/
					alpha.resize(T, vec_dbl(NUM_STATES, 0));
					beta.resize(T, vec_dbl(NUM_STATES, 0));
					gamma.resize(T, vec_dbl(NUM_STATES, 0));
					delta.resize(T, vec_dbl(NUM_STATES, 0));
					psi.resize(T, vec_int(NUM_STATES, 0));
					q.resize(T);
					for (int i = 0; i < NUM_STATES; i++)
					{
						alpha[0][i] = pi[i] * b[i][o[0]];
					}
					for (int t = 1; t < T; t++)
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
					for (int t = T - 2; t >= 0; t--)
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

					int maxp;
					double maxd, tempd;
					for (int i = 0; i < NUM_STATES; i++)
					{
						delta[0][i] = pi[i] * b[i][o[0]];
						psi[0][i] = 0;
					}
					for (int t = 1; t < T; t++)
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
					for (int t = T - 2; t >= 0; t--)
						q[t] = psi[t + 1][q[t + 1]];
					cout << "\npstar: " << pstar << std::endl;
					/*cout << "qstar: ";
					for (auto iq : q)
						cout << iq << ' ';*/
					if ((pstar - pstar_prev) / pstar < 1.0e-4)
					{
						cout << "Converged in " << iter + 1;
						break;
					}
					pstar_prev = pstar;
					std::vector<matrix_dbl> xi(T, matrix_dbl(NUM_STATES, vec_dbl(NUM_STATES, 0)));

					double tot;

					for (int t = 0; t < T - 1; t++)
					{
						tot = 0;
						for (int i = 0; i < NUM_STATES; i++)
						{
							for (int j = 0; j < NUM_STATES; j++)
							{
								xi[t][i][j] = (alpha[t][i] * a[i][j] * b[j][o[t + 1]] * beta[t + 1][j]);
								tot += xi[t][i][j];
							}
						}
						for (int i = 0; i < NUM_STATES; i++)
						{
							for (int j = 0; j < NUM_STATES; j++)
								xi[t][i][j] /= tot;
						}
					}

					for (int t = 0; t < T; t++)
						for (int t = 0; t < T; t++)
						{
							sum = 0;
							for (int i = 0; i < NUM_STATES; i++)
							{
								gamma[t][i] = alpha[t][i] * beta[t][i];
								sum += gamma[t][i];
							}
							for (int i = 0; i < NUM_STATES; i++)
								gamma[t][i] /= sum;
						}

					for (int i = 0; i < NUM_STATES; i++)
					{
						pi[i] = gamma[0][i];
					}
					double trans;
					for (int i = 0; i < NUM_STATES; i++)
					{
						trans = 0;
						for (int t = 0; t < T - 1; t++)
						{
							trans += gamma[t][i];
						}
						for (int j = 0; j < NUM_STATES; j++)
						{
							a[i][j] = 0;
							for (int t = 0; t < T; t++)
							{
								a[i][j] += xi[t][i][j];
							}
							a[i][j] /= trans;
						}
						for (int k = 0; k < M; k++)
						{
							b[i][k] = 0;
							for (int t = 0; t < T; t++)
							{
								if (o[t] == k)
									b[i][k] += gamma[t][i];
							}
							b[i][k] /= trans;
						}
					}
				}
				DEBUGM(a);
				DEBUGM(b);
				DEBUGV(pi);
				DEBUGV(o);
				for (int i = 0; i < NUM_STATES; i++)
				{
					for (int j = 0; j < NUM_STATES; j++)
					{
						if (isnan(a[i][j])) throw std::exception("nan");
						a_avg[i][j] += a[i][j];
					}
				}
				for (int i = 0; i < NUM_STATES; i++)
				{
					for (int j = 0; j < M; j++)
					{
						if (isnan(b[i][j])) throw std::exception("nan");
						b_avg[i][j] += b[i][j];
					}
				}
				for (int i = 0; i < NUM_STATES; i++)
				{
					pi_avg[i] += pi[i];
				}
			}
			for (int i = 0; i < NUM_STATES; i++)
			{
				for (int j = 0; j < NUM_STATES; j++)
				{
					a_avg[i][j] /= 6;
				}
			}
			for (int i = 0; i < NUM_STATES; i++)
			{
				for (int j = 0; j < M; j++)
				{
					b_avg[i][j] /= 6;
				}
			}
			for (int i = 0; i < NUM_STATES; i++)
			{
				pi_avg[i] /= 6;
			}
		}
		DEBUGM(a_avg);
		DEBUGM(b_avg);
		DEBUGV(pi_avg);
		write_output(a_avg, (string)outpath + "a_" + (char)('0' + d));
		write_output(b_avg, (string)outpath + "b_" + (char)('0' + d));
		write_output(pi_avg, (string)outpath + "pi_" + (char)('0' + d));
	}
	return 0;
}

void test()
{
	vec_dbl pi[10];//, pi_init, pi_avg, pi_avg_prev;
	matrix_dbl a[10], b[10], alpha, beta, gamma, delta;//, a_init, b_init, a_avg, b_avg, a_avg_prev, b_avg_prev;
	matrix_int psi;
	vec_int o, q, LBGresult;
	double cur_prob, max_prob, accuracy = 0;
	int pred_dig, total_pred = 0;
	for (int d = 0; d < 10; d++)
	{
		read_mat(a[d], (string)inpath + "a_" + (char)('0' + d), NUM_STATES);
		read_mat(b[d], (string)inpath + "b_" + (char)('0' + d), M);
		read_pi(pi[d], (string)inpath + "pi_" + (char)('0' + d));
	}
	read_obs(LBGresult, (string)inpath + "observations");
	int startpos, endpos;
	for (int test_d = 0; test_d < 10; test_d++)
	{
		for (int file_num = 6; file_num < 10; file_num++)
		{
			max_prob = 0;
			startpos = (test_d*num_of_rec + file_num)*no_of_itern;
			endpos = startpos + no_of_itern;
			o = vec_int(LBGresult.begin() + startpos, LBGresult.begin() + endpos);
			cout << "Digit " << test_d << ", File " << file_num << ":";
			for (int d = 0; d < 10; d++)
			{
				cur_prob = prob(o, a[d], b[d], pi[d]);
				cout<<"\nProb("<<d<<")="<<cur_prob;
				if (cur_prob > max_prob)
				{
					max_prob = cur_prob;
					pred_dig = d;
				}
			}
			cout << " Predicted " << pred_dig << '\n';
			if (pred_dig == test_d)
				accuracy += 1;
			total_pred++;
		}
	}
	cout << "Accuracy: " << 100 * accuracy / total_pred;
}

int main()
{
	train();
	test();
	return 0;
}
