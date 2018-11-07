
// < Change these as required
constexpr auto inpath = "./";						//Path of the recorded files
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
constexpr auto no_of_itern = 5;					//Number of windows in each recording
constexpr auto K = 5;							//Number of K-Means cluster
constexpr auto convergence_r = 1.0e-20;			//Convergence radius
constexpr int weight[12] = { 1,3,7,13,19,22,25,33,42,50,56,61 };
