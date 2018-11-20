#include "pch.h"

vec_dbl hamming(vec_dbl s, int wsize)
{
	if (wsize < 0)
		wsize = s.size();
	for (int i = 0; i < wsize; i++)
	{
		s[i] = s[i] * (0.56 - 0.46*cos((2 * M_PI * i) / wsize));
	}
	return s;
}

vec_dbl calc_R(vec_dbl s, int p)
{
	vec_dbl R(p + 1, 0);
	for (int i = 0; i <= p; i++)
	{
		for (int m = 0; i + m < s.size(); m++)
			R[i] += s[m] * s[i + m];
	}
	return R;
}

vec_dbl calc_A(vec_dbl &R)
{
	int p = R.size();
	double k;
	vec_dbl A(p, 0), E(p, 0), A_temp;
	E[0] = R[0];

	for (int i = 1; i < p; i++)
	{
		A_temp = A;
		k = R[i];
		for (int j = 1; j < i; j++)
			k -= A_temp[j] * R[i - j];
		k /= E[i - 1];
		A[i] = k;
		for (int j = 1; j < i; j++)
			A[j] = A_temp[j] - k * A_temp[i - j];
		E[i] = (1 - k * k)*E[i - 1];
	}
	return A;
}

vec_dbl calc_C(vec_dbl &A, double sigma_sq)
{
	int p = A.size();
	vec_dbl C(A);
	C[0] = log2(sigma_sq);
	for (int i = 1; i < p; i++)
	{
		for (int m = 1; m < i; m++)
			C[i] += m * C[m] * A[i - m] / i;
	}
	return C;
}

double dist(vec_dbl ref, vec_dbl target)
{
	double res = 0, diff;
	for (int i = 0; i < ref.size(); i++)
	{
		diff = ref[i] - target[i];
		res += weight[i] * diff*diff;
	}
	return res;
}

vec_dbl readFile(string inputfile)
{
	std::ifstream fin;
	vec_dbl v(0);

	fin.open(inputfile);

	if (fin.fail())
	{
		std::cout << "Unable to open \"" << inputfile << "\"!\n\n";
		std::cin.ignore(1000, '\n');
		return v;
	}

	string tempstr;

	while (true)				//SKIP INITIAL LINES
	{
		std::getline(fin, tempstr);
		if (tempstr[0] == '-' || tempstr[0] < 58 && tempstr[0]>47)
			break;
	}

	double inp = stod(tempstr);
	int max = 1;
	while (!fin.eof())			//READ SAMPLES
	{
		inp -= DC;				//DC shift
		v.push_back(inp);
		if (abs(inp) > max)
			max = inp;
		fin >> inp;
	}
	double nrm = (double)5000 / max;	//Normalisation factor
	for (int i = 0; i < v.size(); i++)
	{
		v[i] *= nrm;
	}
	fin.close();
	return v;
}

std::pair<int, int> VAD(vec_dbl v)
{
	double e = 0;
	int i;
	for (i = 0; i < 400; i++)	//CALCULATE ENERGY RANGE FOR SILENCE(NOISE)
	{
		e += v[i] * v[i];
	}
	double et = 4 * e;			//SET ENERGY THRESHHOLD

	for (; i < v.size(); i++)	//FIND START POINT
	{
		e += v[i] * v[i] - v[i - 400] * v[i - 400];
		if (e > et)
			break;
	}
	int start = i - 400;
	int end = v.size();
	//for (; i < v.size(); i++)	//FIND END POINT
	//{
	//	e += v[i] * v[i] - v[i - 400] * v[i - 400];
	//	if (e < et)
	//	{
	//		end = i - 400;
	//		break;
	//	}
	//}
	////std::cout << start << " to " << end << " , total frames = "<< end-start << '\n';
	return std::make_pair(start, end);
}

matrix_dbl calc_C_for_file(string inputfile)
{
	matrix_dbl C(0);
	vec_dbl R, A, v = readFile(inputfile);
	if (v.size() < 1)
		return C;

	C.resize(no_of_itern);

	int start = 0, end = v.size();
	std::pair<int, int> tempp = VAD(v);
	start = tempp.first; end = tempp.second;
	//start += (end - start) / 2 - (window_size + stride * (no_of_itern - 1)) / 3;

	for (int itern = 0; itern < no_of_itern; itern++)
	{
		vec_dbl samples(v.begin() + start + stride * itern, v.begin() + start + stride * itern + window_size);

		//R = calc_R(samples, p);
		R = calc_R(hamming(samples), p);

		A = calc_A(R);

		C[itern] = calc_C(A, R[0]);
	}
	return C;
}

void makeUniverse(string outfile)
{
	string infile;
	std::ofstream fout(outfile);
	char dig[] = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' };

	matrix_dbl refC[10], C;
	for (int i = 0; i < 10; i++)
	{
		for (int j = 0; j < num_of_rec; j++)
		{
			infile = inpath + (string)"150101012_" + dig[i] + (string)"_" + std::to_string(j + 1) + (string)".txt";
			//std::cout << inputfile << " : ";
			C = calc_C_for_file(infile);
			refC[i].insert(refC[i].end(), C.begin(), C.end());
			for (int i = 0; i < no_of_itern; i++)
			{
				fout << C[i][1];
				for (int j = 2; j <= p; j++)
					fout << ',' << C[i][j];
				fout << '\n';
			}
		}
	}
	//return C;
}

matrix_dbl readUniverse(string infile)
{
	matrix_dbl C(0);
	vec_dbl C_inp(p);
	std::ifstream fin(infile);
	char comma_in;
	while (true)
	{
		fin >> C_inp[0];
		if (fin.eof())
			break;
		for (int i = 1; i < p; i++)
			fin >> comma_in >> C_inp[i];
		C.push_back(C_inp);
	}
	return C;
}

matrix_dbl calc_means(const matrix_dbl &C, const vec_int &class_, int K)
{
	matrix_dbl means(K, vec_dbl(p, 0));
	vec_int count(K, 0);
	for (int i = 0; i < class_.size(); i++)
	{
		for (int j = 0; j < p; j++)
			means[class_[i]][j] += C[i][j];
		count[class_[i]]++;
	}
	for (int i = 0; i < K; i++)
		for (int j = 0; j < p; j++)
			if (count[i]) means[i][j] /= count[i];
	return means;
}

void classify(const matrix_dbl &C, const matrix_dbl &means, vec_int &class_, int K)
{
	double d;
	for (int i = 0; i < C.size(); i++)
	{
		int jmin = -1;
		double min = std::numeric_limits<double>::max();
		for (int j = 0; j < K; j++)
		{
			d = dist(C[i], means[j]);
			if (d < min)
			{
				min = d;
				jmin = j;
			}
		}
		class_[i] = jmin;
	}
}

void update_means(const matrix_dbl &C, matrix_dbl &means, vec_int &class_, const int K)
{
	matrix_dbl prev_means = means;
	double diff;
	bool break_loop;
	while (true)
	{
		break_loop = true;
		classify(C, means, class_, K);
		means = calc_means(C, class_, K);
		for (int i = 0; i < K; i++)
		{
			diff = dist(means[i], prev_means[i]);
			std::cout << diff << ' ';
			if (diff > convergence_r)
			{
				break_loop = false;
			}
		}
		std::cout << '\n';
		prev_means = means;
		if (break_loop)
			break;
	}
}

vec_int kmeansClassifier(const matrix_dbl &C, int K)
{
	int n = C.size();
	vec_int class_(n);
	matrix_dbl means(K, vec_dbl(p));
	srand(time(NULL));
	for (int i = 0; i < n; i++)
		class_[i] = rand() % K;

	means = calc_means(C, class_, K);
	/*double min, max;
	for (int i = 0; i < p; i++)
	{
		min = C[std::min_element(C.begin(), C.end(), [i](const vec_dbl a, const vec_dbl b) {return a[i] < b[i]; }) - C.begin()][i];
		max = C[std::max_element(C.begin(), C.end(), [i](const vec_dbl a, const vec_dbl b) {return a[i] < b[i]; }) - C.begin()][i];
		for (int j = 0; j < K; j++)
		{
			means[j][i] = min + ((double)rand() / RAND_MAX)*(max - min);
		}
	}*/
	update_means(C, means, class_, K);
	classify(C, means, class_, K);
	std::ofstream fout(kMeansCodebookFile);
	for (int i = 0; i < K; i++)
	{
		fout << means[i][0];
		for (int j = 1; j < p; j++)
			fout << ',' << means[i][j];
		fout << '\n';
	}
	return class_;
}

vec_dbl calc_stddev(matrix_dbl C, vec_dbl mean)
{
	vec_dbl stddev(p);
	for (int i = 0; i < p; i++)
	{
		for (int j = 0; j < C.size(); j++)
			stddev[i] += (C[j][i] - mean[i])*(C[j][i] - mean[i]);
		stddev[i] /= C.size();
		stddev[i] = sqrt(stddev[i]);
	}
	return stddev;
}

vec_int LBGClassifier(matrix_dbl C)
{
	matrix_dbl means, means_copy;
	vec_int class_(C.size(), 0);
	means.push_back(calc_means(C, class_, 1)[0]);
	vec_dbl epsilon = calc_stddev(C, means[0]);
	for (auto &x : epsilon)
		x /= 1000;
	while (means.size() < K)
	{
		//means=split(means,epsilon);
		means_copy = means;
		for (auto &x : means)
			for (int i = 0; i < p; i++)
				x[i] += epsilon[i];
		for (auto &x : means_copy)
			for (int i = 0; i < p; i++)
				x[i] -= epsilon[i];
		means.insert(means.end(), means_copy.begin(), means_copy.end());
		update_means(C, means, class_, means.size());
	}
	classify(C, means, class_, means.size());
	std::ofstream fout(LBGCodebookFile);
	for (int i = 0; i < means.size(); i++)
	{
		fout << means[i][0];
		for (int j = 1; j < p; j++)
			fout << ',' << means[i][j];
		fout << '\n';
	}
	return class_;
}

// In general, ignore this file, but keep it around if you are using pre-compiled headers.
