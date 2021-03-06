#include "Utils.hpp"
#include "Configurator.hpp"
#include <cstring>
#include <fstream>
#include <NTL/ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/LLL.h>
#include <random>
#include <cstdlib>
#include <algorithm>

using namespace NTL;
using namespace std;

Pruning to_prune(const std::string& s);
int readConfig(string path);

void searchShortestSeedCandidates(int dim, int startseed, int endseed) {
	double Aest[11];
	int seedest[11];

	for(int i=0; i < 11; i++) {
		Aest[i] = numeric_limits<double>::max();
		seedest[i] = -1;
	}

	for(int seed = startseed; seed <= endseed; seed++) {
		std::string filename = "svpc-e1-d";
		filename += std::to_string(dim);
		filename += "-";
		filename += std::to_string(seed);
		filename += ".txt";

		std::string command = "/home/mb31simi/svpgenerator/generate_random --dim ";
		command += std::to_string(dim);
		command += " --seed ";
		command += std::to_string(seed);
		command += " --bit ";
		command += std::to_string(dim*10);
		command += " --bitexponent 1.0 &> /home/mb31simi/NewPenum/penum/tinput/";
		command += filename;

		cout << command << endl;
		std::system(command.c_str());

		// Read Matrix
	    NTL::mat_ZZ B;
	    std::ifstream basestream;
	    std::string filepath = "tinput/" + filename;
	    std::cout << "Trying to read: " << filepath << std::endl;
	    basestream.open(filepath.c_str(), ios_base::in);
	    try {
	    	basestream >> B;
	    }
	    catch (std::exception & e) {
	    	cout << "Error while reading basis. Exiting" << std::endl;
	    	exit(-1);
	    }
	    basestream.close();

	    std::cout << "Read a basis of dimension: "
	    		<< B.NumRows() << " / " << B.NumCols()
				<< " with |b_0| = " << lengthOfVector<ZZ, RR> (B[0])
				<< std::endl;

	    LLL_XD(B, 0.6);

		double gauss = calcGaussHeuristic<ZZ>(B);
			cout << "Gaussian heuristic for lambda_1(b): "
					<< sqrt(gauss) << endl;

		// Check if it is a candidate for the hitlist
		int pos = -1;
		for(int i=0; i<10; i++) {
			if(gauss < Aest[i]) {
				pos = i;
				cout << "New min at position " << i << endl;
				break;
			}
		}

		// Adapt the value and seed array
		if(pos > -1) {
			for(int i=9; i>pos; i--) {
				//cout << "Moving " << i <<
				Aest[i] = Aest[i-1];
				seedest[i] = seedest[i-1];
			}
			Aest[pos] = gauss;
			seedest[pos] = seed;

			// Show the new situation
			cout << "Candidatelist changed:" << endl;

			for(int i=0; i<10;i++) {
				cout << "Pos" << i+1 << ": " << Aest[i] << "[seed=" << seedest[i] << "]" << endl;
			}
		}
	}
}

Pruning to_prune(const std::string& s)
{
    if (s == "linear" || s == "Linear" || s == "LINEAR") return LINEAR;
    if (s == "NO" || s == "No" || s == "no") return NO;
    if (s == "soft") return SOFT;
    if (s == "STEP" || s == "Step" || s == "step") return STEP;
    if (s == "PARLINEAR" || s == "Parlinear" || s == "parlinear") return PARLINEAR;
    if (s == "EXTREME" || s == "Extreme" || s == "extreme") return EXTREME;
    throw std::runtime_error("invalid conversion from text to periods");
}


int readPruning(string path) {
	std::ifstream file(path.c_str());

	if (!file)
	{
	    cerr << "No valid pruning function file at " << path << endl;
	    return -1;
	}

	std::string   line;
	double dval;
	int valcnt = 0;

	while(std::getline(file, line)) {
		valcnt ++;
	}	
	file.close();	
	cout << "Reading a pruning function of dimension " << valcnt << endl;
	Configurator::getInstance().ext_pruning_function.clear();
	Configurator::getInstance().ext_pruning_function.resize(valcnt);
	Configurator::getInstance().ext_pruning_function.reserve(valcnt);

	file.open(path.c_str());
	int pos = valcnt-1;
	while(std::getline(file, line)) {
		std::stringstream   linestream(line);
		linestream >> dval;
		//cout << "dval[" << pos << "] = " << dval << endl;
		Configurator::getInstance().ext_pruning_function[pos] = dval;
		pos--;
	}

	return 0;
}

int readConfig(string path) {
	std::string locpath = path;
	locpath += "penum.conf";

	std::ifstream file(locpath.c_str());

	if (!file)
	{
	    cerr << "No valid conf file at " << locpath << endl;
	    return -1;
	}

	std::string   line;

	while(std::getline(file, line))
	{
		std::stringstream   linestream(line);
		std::string         data;

		string command;
		string str_value;
		int int_value;
		double double_value;

		linestream >> command;
		if(command.length()>0) {
			// Skip comments
			if(command[0] == '#')
				continue;
		}

		// Parse commands an register in configurator
		// Should the final enumeration call be pruned, if yes, which pruning will be activated?
        if(command == "--enumpruning") {
        	try {
        		linestream >> str_value;
        		Configurator::getInstance().enum_prune = to_prune(str_value);
        		cout << "Activating " << str_value << " pruning." << endl;
        	}

        	catch (exception& e) {
        		std::cerr << "Missing enumpruning in argument. Terminating." << std::endl;
        	}

        }

        if(command == "--pruneparam") {
        	try {
        		linestream >> double_value;
        		Configurator::getInstance().prune_param = double_value;
        		cout << "Setting pruning parameter to " << double_value << "." << endl;
        	}

        	catch (exception& e) {
        		std::cerr << "Missing pruning parameter in argument." << std::endl;
        	}

        }


        if(command == "--candsize") {
        	try {
        		linestream >> int_value;
        		Configurator::getInstance().cand_queue_size = int_value;
        		cout << "Setting candidate-queue size to " << int_value << "." << endl;
        	}

        	catch (exception& e) {
        		std::cerr << "Missing candsize in argument." << std::endl;
        	}

        }

        if(command == "--bkzheight") {
        	try {
        		linestream >> int_value;
        		Configurator::getInstance().bkz_height = int_value;
        		cout << "Setting bkz-height to " << int_value << "." << endl;
        	}

        	catch (exception& e) {
        		std::cerr << "Missing bkz-height in argument." << std::endl;
        	}

        }

        if(command == "--parthres") {
        	try {
        		linestream >> int_value;
        		Configurator::getInstance().par_threshold = int_value;
        		cout << "Setting parallel-enum threshold to dimension " << int_value << "." << endl;
        	}

        	catch (exception& e) {
        		std::cerr << "Missing parthres in argument. Terminating." << std::endl;
        	}

        }
        
        if(command == "--BKZinstances") {
        	try {
        		linestream >> int_value;
        		Configurator::getInstance().BKZinstances = int_value;
        		cout << "Setting parallel-BKZ instances to " << int_value << "." << endl;
        	}

        	catch (exception& e) {
        		std::cerr << "Missing BKZinstances in argument. Terminating." << std::endl;
        	}
        }
        
        if(command == "--prebeta") {
        	try {
        		linestream >> int_value;
        		Configurator::getInstance().prebeta = int_value;
        		cout << "Setting pre-beta to " << int_value << "." << endl;
        	}

        	catch (exception& e) {
        		std::cerr << "Missing pre-beta in argument. Terminating." << std::endl;
        	}
        }

        if(command == "--ntlbkz") {
        	try {
        		Configurator::getInstance().use_fplll_bkz = false;
        		cout << "Enable usage of NTL BKZ" << endl;
        	}

        	catch (exception& e) {
        		std::cerr << "Error concerning setting of NTL BKZ." << std::endl;
        	}

        }
        

        if(command == "--trials") {
        	try {
        		linestream >> int_value;
        		Configurator::getInstance().trials = int_value;
        		cout << "Setting randomized trials to " << int_value << "." << endl;
        	}

        	catch (exception& e) {
        		std::cerr << "Missing randomize trial count in argument." << std::endl;
        	}

        }

        if(command == "--nogaussestimate") {
        	try {
        		Configurator::getInstance().do_gauss_estimate = false;
        		cout << "Disabling usage of Gauss heuristic" << endl;
        	}

        	catch (exception& e) {
        		std::cerr << "Error concerning setting of Gauss heuristic." << std::endl;
        	}

        }

        if(command == "--bkzautoheight") {
        	try {
        		Configurator::getInstance().do_bkz_auto_tuning = true;
        		cout << "Enabling BKZ-auto-height adaption." << endl;
        	}

        	catch (exception& e) {
        		std::cerr << "Error concerning setting of BKZ-auto-height adaption." << std::endl;
        	}

        }




		// Parse input arguments
	    /*for(int i = 1; i < argc; i++) {
	    	std::string input = std::string(argv[i]);
	        if(input == "--vecfile") {
	            if(argc <= i) {
	                std::cerr << "Missing vectorfile in argument. Terminating." << std::endl;
	                exit(-1);
	            }
	            cout << "Only executing vector length mode." << endl;
	            vecpath = std::string(argv[i+1]);
	            i++;
	        }
	        if(input == "--beta") {
	            if(argc <= i) {
	                std::cerr << "Missing beta in argument. Terminating." << std::endl;
	                exit(-1);
	            }
	            beta = atoi(argv[i+1]);
	            Configurator::getInstance().glob_beta = beta;
	            i++;
	        }
	        if(input == "--sheight") {
	            if(argc <= i+1) {
	                std::cerr << "Missing sheight in argument. Terminating." << std::endl;
	                exit(-1);
	            }
	            candidate_height = atoi(argv[i+1]);
	            Configurator::getInstance().forced_height = candidate_height;
	            i++;
	        }
	        if(input == "--noautoheight") {
	        	Configurator::getInstance().auto_height_tuning = false;
	        }
	        if(input == "--delta") {
	            if(argc <= i) {
	                std::cerr << "Missing delta in argument. Terminating." << std::endl;
	                exit(-1);
	            }
	            delta = atof(argv[i+1]);
	            Configurator::getInstance().glob_delta = delta;
	            i++;
	        }
	        if(input == "--mbenum") {
	            Configurator::getInstance().use_ntl_enum = false;
	        }
	        if(input == "--onlyLLLexport") {
	            if(argc <= i) {
	                std::cerr << "Missing output file argument for LLL. Terminating." << std::endl;
	                exit(-1);
	            }
	            reducedfilepath = std::string(argv[i+1]);
	            mode = 2;
	            i++;
	        }
	        if(input == "--noLLL") {
	        	dolll = false;
	        }
	    }*/


	}

	cout << "Read conf file at " << locpath << endl;
	return 0;
}
int intRandTS(const int & min, const int & max) {
    static thread_local std::mt19937 generator;
    std::uniform_int_distribution<int> distribution(min,max);
    return distribution(generator);
}


mat_ZZ randomizeMatrix(const mat_ZZ A, const vec_ZZ randoms, const int dim, const int seed){
       ZZ one;
       one = 1;
       mat_ZZ newmat, P;
       newmat.SetDims(dim,dim);
       P.SetDims(dim,dim);
       int newmati, newmatj, randompicker;

       srand(seed * time(NULL));

       for (int i=0; i<((dim*dim-dim)/2); i++){
               //newmati = rand()%(dim-1);
               //newmatj = rand()%(newmati+1);
               //randompicker = rand()%5;

    	   newmati = intRandTS(0, dim-2);
    	   newmatj = intRandTS(0, newmati);
    	   randompicker = intRandTS(0, 4);

               newmat[newmati+1][newmatj] = randoms[randompicker];
       }
       for (int i=0; i<dim; i++){
               newmat[i][i] = one;
               P[i][i] = one;
       }
       mul(newmat, newmat, A);
       int first,second;
       ZZ temp;
       for(int i=0; i<dim;i++){
               first = rand()%dim;
               second = rand()%dim;
               for (int j=0; j<dim; j++){
                       temp = P[j][first];
                       P[j][first] = P[j][second];
                       P[j][second] = temp;
               }
       }
       mul(newmat, P, newmat);
       return newmat;
}

mat_RR GSO(const mat_ZZ& input, mat_RR& mu, vec_RR& c) {
	ComputeGS(input, mu, c);
	// set the diagonal
	for (int i = 0; i < input.NumRows(); ++i)
		mu[i][i] = RR(1);

	mat_RR Bstar = inv(mu) * conv < mat_RR > (input);

	if(Configurator::getInstance().do_gauss_estimate) {
		Configurator::getInstance().gauss_estimate = (calcGaussHeuristic<RR>(Bstar));
		cout << "Gaussian heuristic for lambda_1(B*): "
				<< Configurator::getInstance().gauss_estimate << endl;
	}



	return Bstar;
}

void sortBasesByQuality(const std::vector<mat_ZZ>& Bcands, std::vector<int>& order) {
	int ts = Bcands.size();
	cout << "Sorting " << ts << " bases." << endl;
	vector<double> len0;
	vector<int> num;

	// Just fill with 1 ... dimension
	vector<double> x;
	for(int j=0; j < Bcands[0].NumRows(); j++) {
		x.push_back(double(j));
	}

	len0.reserve(ts);
	len0.resize(ts);

	// Try it with linear regression
	for (int tid=0; tid < ts; tid ++) {
		mat_RR mu;
		vec_RR bstar;
		GSO(Bcands[tid], mu, bstar);

		vector<double> y;

		for(int j=0; j < Bcands[tid].NumRows(); j++) {
			y.push_back(conv<double>(log10(bstar[j])));
		}

		if (x.size() != y.size()) {
			cerr << "[sortBasesByQuality]: Size mismatch of arrays" << endl;
		}
		double n = x.size();

		double avgX = accumulate(x.begin(), x.end(), 0.0) / n;
		double avgY = accumulate(y.begin(), y.end(), 0.0) / n;

		double numerator = 0.0;
		double denominator = 0.0;

		for (int i = 0; i<n; ++i) {
			numerator += (x[i] - avgX) * (y[i] - avgY);
			denominator += (x[i] - avgX) * (x[i] - avgX);
		}

		if (denominator == 0) {
			cerr << "[sortBasesByQuality]: denominator is zero" << endl;
		}
		//cout << numerator / denominator << endl;
		len0[tid] = numerator / denominator;
	}

	// Sort by length and adapt order of processing
	for(int i=0; i < ts; i++) {
		double minval = len0[order[i]];
		int minpos=i;
		for(int j=i; j < ts; j++) {
			if(len0[order[j]] < minval) {
				minval = len0[order[j]];
				minpos = j;
			}
		}
		// Swap min
		int old = order[i];
		order[i] = order[minpos];
		order[minpos] = old;
	}

	std::reverse(order.begin(), order.end());
	for(int i=0; i<ts; i++) {
			cout << len0[order[i]] << " ";
		}
		cout << endl;





	/*for(int i=0; i<ts; i++) {
		len0.clear();
		num.clear();
		len0.reserve(Bcands[i].NumRows());
		len0.resize(Bcands[i].NumRows());
		num.reserve(Bcands[i].NumRows());
		num.resize(Bcands[i].NumRows());

		mat_RR mu;
		vec_RR bstar;
		GSO(Bcands[i], mu, bstar);

		for(int j=0; j < Bcands[i].NumRows(); j++) {
			len0[j] = conv<double>(log(bstar[j]));
			//cout << len0[j] << " ";
			num[j] = j;
		}
		//cout << endl;

		const auto n    = len0.size();
		const auto s_x  = std::accumulate(num.begin(), num.end(), 0.0);
		const auto s_y  = std::accumulate(len0.begin(), len0.end(), 0.0);
		const auto s_xx = std::inner_product(num.begin(), num.end(), num.begin(), 0.0);
		const auto s_xy = std::inner_product(num.begin(), num.end(), len0.begin(), 0.0);
		const auto a    = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
		cout << "Slope is:" << a << endl;

	}
	cout << endl;



	// Get length of each b_0*
	len0.reserve(ts);
	len0.resize(ts);
	num.reserve(ts);
	num.resize(ts);
	for(int i=0; i<ts; i++) {
		len0[i] = conv<double>(lengthOfVector<ZZ, RR> (Bcands[i][0]));
		num[i] = i;
	}

	// Sort by length and adapt order of processing
	for(int i=0; i < ts; i++) {
		double minval = len0[order[i]];
		int minpos=i;
		for(int j=i; j < ts; j++) {
			if(len0[order[j]] < minval) {
				minval = len0[order[j]];
				minpos = j;
			}
		}
		// Swap min
		int old = order[i];
		order[i] = order[minpos];
		order[minpos] = old;
	}*/






	return;
}
