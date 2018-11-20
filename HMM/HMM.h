
// < Change these as required
constexpr auto inpath = "./";						//Path of the recorded files
constexpr auto outpath = "./";						//Path of the recorded files
constexpr auto num_of_rec = 10;						//Number of recordings for each vowel
constexpr auto CreateUniverseFile = true;			//Set false to use existing Universe.csv
constexpr auto UniverseFile = "./Universe.csv";
constexpr auto OutputFile = "./result.csv";
constexpr auto kMeansCodebookFile = "./means_kMeans.csv";
constexpr auto LBGCodebookFile = "./means_LBG.csv";
// >


constexpr auto DC = 0;							//DC Offset (Obtained to be approx. 0 for my setup)
constexpr auto p = 12;							//Number of Cepstral Coeff.
constexpr auto window_size = 320;
constexpr auto stride = 80;						//25% stride
constexpr auto no_of_itern = 40;					//Number of windows in each recording
constexpr auto K = 16;							//Number of K-Means cluster
constexpr auto convergence_r = 1.0e-20;			//Convergence radius
constexpr int weight[12] = { 1,3,7,13,19,22,25,33,42,50,56,61 };

constexpr int NUM_STATES = 5;
constexpr int M = 16;

using vec_int = std::vector<int>;
using vec_dbl = std::vector<double>;
using matrix_dbl = std::vector<vec_dbl>;
using std::string;
using std::cout;
using std::cin;

vec_dbl hamming(vec_dbl s, int wsize = -1);
vec_dbl calc_R(vec_dbl s, int p);
vec_dbl calc_A(vec_dbl &R);
vec_dbl calc_C(vec_dbl &A, double sigma_sq);
double dist(vec_dbl ref, vec_dbl target);
vec_dbl readFile(string inputfile);
std::pair<int, int> VAD(vec_dbl v);
matrix_dbl calc_C_for_file(string inputfile);
void makeUniverse(string outfile = UniverseFile);
matrix_dbl readUniverse(string infile = UniverseFile);
matrix_dbl calc_means(const matrix_dbl &C, const vec_int &class_, int K = ::K);
void classify(const matrix_dbl &C, const matrix_dbl &means, vec_int &class_, int K = ::K);
void update_means(const matrix_dbl &C, matrix_dbl &means, vec_int &class_, const int K = ::K);
vec_int kmeansClassifier(const matrix_dbl &C, int K = ::K);
vec_dbl calc_stddev(matrix_dbl C, vec_dbl mean);
vec_int LBGClassifier(matrix_dbl C);