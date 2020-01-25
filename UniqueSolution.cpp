
// UniqueSolution.cpp : The Number of Unique Identities Created by the Six Major Operations on Any Two Integers 0-X,
//     Function Finder
//
#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <cmath>
#include<math.h>
#include <vector>


std::vector<std::vector<std::vector<double>>> uniqueValues(6);


void calcSeq(double depth) {
	for (double it = 0; it < 6; it++) {
		uniqueValues[it].resize(depth);
	}
	for (double i = 0; i < 6; i++) {
		for (double it = 0; it < depth; it++) {
			for (double iu = 0; iu <= it; iu++) {
				for (double ii = 0; ii <= it; ii++) {
					if (i == 0) {
						bool found=false;
						for (int io = 0; io < uniqueValues[i][it].size(); io++) {
							if (!ii == 0) {
								if (uniqueValues[i][it][io] == std::pow(iu, 1.0 / ii)) {
									found=true;
								}
							}
							if (ii == 0) {
								if (iu!=0 && uniqueValues[i][it][io] ==(M_PI)) {
									found=true;
								}
								else {
									if (iu == 0 && uniqueValues[i][it][io] == 1) {
										found = true;
									}
								}
							}
						}
						if (!found) {
							
							if (!ii == 0) {
								uniqueValues[i][it].push_back(std::pow(iu, 1.0 / ii));
							}
							if (ii == 0) {
								if (iu == 0) {
									uniqueValues[i][it].push_back(1.0);
								}
								else {
									uniqueValues[i][it].push_back(M_PI);
								}
							}
						}
					}
					if (i == 1) {
						bool found = false;
						for (int io = 0; io < uniqueValues[i][it].size(); io++) {
							if (!ii == 0) {
								if (uniqueValues[i][it][io] == (iu / ii)) {
									found=true;
								}
							}
							if (ii == 0) {
								if (uniqueValues[i][it][io] == (DBL_MAX)) {
									found=true;
								}
							}
						}
						if (!found) {

							if (!ii == 0) {
								uniqueValues[i][it].push_back(iu/ii);
							}
							if (ii == 0) {
								uniqueValues[i][it].push_back((DBL_MAX));
							}
						}
					}
					if (i == 2) {
						bool found = false;
						for (int io = 0; io < uniqueValues[i][it].size(); io++) {
							if (uniqueValues[i][it][io] == (iu-ii)) {
								found=true;
							}
						}
						if (!found) {
							uniqueValues[i][it].push_back(iu - ii);
						}
					}
					if (i == 3) {
						bool found = false;
						for (int io = 0; io < uniqueValues[i][it].size(); io++) {
							if (uniqueValues[i][it][io] == (iu + ii)) {
								found=true;
							}
						}
						if (!found) {
							uniqueValues[i][it].push_back(iu + ii);
						}
					}
					if (i == 4) {
						bool found = false;
						for (int io = 0; io < uniqueValues[i][it].size(); io++) {
							if (uniqueValues[i][it][io] == (iu * ii)) {
								found=true;
							}
						}
						if (!found) {
							uniqueValues[i][it].push_back(iu * ii);
						}
					}
					if (i == 5) {
						bool found = false;
						for (int io = 0; io < uniqueValues[i][it].size(); io++) {
							if (uniqueValues[i][it][io] == std::pow(iu, ii)) {
								found=true;
							}
						}
						if (!found) {
							uniqueValues[i][it].push_back(std::pow(iu, ii));
						}
					}
				}//root,divide,subtract,add,multiply,pow
			}//root,divide,subtract,add,multiply,pow
		}//root,divide,subtract,add,multiply,pow
		std::cout << ".";
	}//root,divide,subtract,add,multiply,pow
	std::cout << "\n \n \n Done. \n Roots sequence: \n";
	for (int i= 0; i < depth; i++) {
		std::cout << uniqueValues[0][i].size()<<", ";
	}
	std::cout << "\n Quotient sequence: \n";
	for (int i = 0; i < depth; i++) {
		std::cout << uniqueValues[1][i].size() << ", ";
	}
	std::cout << "\n Difference sequence: \n";
	for (int i = 0; i < depth; i++) {
		std::cout << uniqueValues[2][i].size() << ", ";
	}
	std::cout << "\n Addend sequence: \n";
	for (int i = 0; i < depth; i++) {
		std::cout << uniqueValues[3][i].size() << ", ";
	}
	std::cout << "\n Product sequence: \n";
	for (int i = 0; i < depth; i++) {
		std::cout << uniqueValues[4][i].size() << ", ";
	}
	std::cout << "\n Powers sequence: \n";
	for (int i = 0; i < depth; i++) {
		std::cout << uniqueValues[5][i].size()<<", ";
	}//root,divide,subtract,add,multiply,pow
}//root,divide,subtract,add,multiply,pow






bool checkCloser(std::vector<double>graph, std::vector<double>pregraph, std::vector<double>sequence, double factorDepth) {
	std::vector<std::vector<bool>>nearingAnswerVector(2);
	nearingAnswerVector.resize(0);
	nearingAnswerVector.resize(2);
	nearingAnswerVector[0].resize(0);
	nearingAnswerVector[1].resize(0);
	for (int ic = 0; ic < sequence.size(); ic++) {
		if (abs(graph.at(sequence[ic]*factorDepth))<= abs(pregraph.at(sequence[ic] * factorDepth))) {
			for (int f = 0; f < factorDepth; f++) {
				nearingAnswerVector[0].push_back(true);
			}
		}
		if (abs(graph.at(sequence[ic] * factorDepth)) > abs(pregraph.at(sequence[ic] * factorDepth))) {
			for (int f = 0; f < factorDepth; f++) {
				nearingAnswerVector[1].push_back(true);
			}
		}
	}
	
	for (int ib = 0; ib < (sequence.size()) * factorDepth; ib++) {//maybe seqsize-1
		if (abs(pregraph[ib]) < abs(graph[ib])) {
			nearingAnswerVector[0].push_back(true);
		}
		if (abs(graph[ib]) <= abs(pregraph[ib])) {
			nearingAnswerVector[1].push_back(true);
		}
		for (int in = 0; in < sequence.size();in++) {
			if (ib / factorDepth == sequence[in]) {
				if (abs(pregraph[ib]) > abs(graph[ib])) {
					nearingAnswerVector[0].resize(nearingAnswerVector[0].size()-1);
				}
				if (abs(graph[ib]) <= abs(pregraph[ib])) {
					nearingAnswerVector[1].resize(nearingAnswerVector[1].size() - 1);
				}
			}
		}
	}
	if (nearingAnswerVector[0].size() > nearingAnswerVector[1].size()) {
		//std::cout << "\n found a factor...";
		return true;
	}
	if (nearingAnswerVector[0].size() <= nearingAnswerVector[1].size()) {
		//std::cout << "\n did not find a factor...";
		return false;
	}
	
}



void modulaSeq(std::vector<std::vector<std::vector<double>>>seq, double factorDepth) {
	double depth;
	std::vector<double>sequence;
	std::vector<double>graph;
	std::vector<double>pregraph;
	std::vector<std::vector<double>>frequencies;
	
	
	frequencies.resize(0);
	frequencies.resize(6);

	
		

	std::cout << "\n Searching for combinations of cosine modulations to fit sequence, scanning range six times for patterning.";
	for (int ik = 0; ik < 6; ik++) {
		std::cout << "\n " << ik << "Set processing cosines, \n";
		sequence.resize(0);
		depth = seq[ik][seq[ik].size()-1].size();
		for (int it = 0; it < seq[ik].size(); it++) {
			sequence.push_back(seq[ik][it].size());
		}
		graph.resize(0);
		pregraph.resize(0);
		graph.resize(factorDepth * (depth+1));
		pregraph.resize(factorDepth * (depth+1));
		for (int it = 0; it < factorDepth; it++) {
			for (int i = 0; i < graph.size(); i++) {
				pregraph[i] += cos((M_PI *.5* i  * it / (factorDepth)));
			}
			if (checkCloser(graph, pregraph, sequence, factorDepth)) {
				for (int i = 0; i < graph.size(); i++) {
					graph[i] += cos((M_PI * i * .5 * it / (factorDepth)));
				}
				frequencies[ik].push_back(it / (factorDepth));
			}
			if (!checkCloser(graph, pregraph, sequence, factorDepth)) {
				for (int i = 0; i < graph.size(); i++) {
					pregraph[i] -= cos((M_PI * i * .5 * it / (factorDepth)));
				}
			}
		}
		for (int it = 0; it < factorDepth; it++) {
			for (int i = 0; i < graph.size(); i++) {
				pregraph[i] += cos((M_PI * i * .5 * it / (factorDepth)));
			}
			if (checkCloser(graph, pregraph, sequence, factorDepth)) {
				for (int i = 0; i < graph.size(); i++) {
					graph[i] += cos((M_PI * i * .5 * it / (factorDepth)));
				}
				frequencies[ik].push_back(it / (factorDepth));
			}
			if (!checkCloser(graph, pregraph, sequence, factorDepth)) {
				for (int i = 0; i < graph.size(); i++) {
					pregraph[i] -= cos((M_PI *.5 *i * it / (factorDepth)));
				}
			}
		}
		for (int it = 0; it < factorDepth; it++) {
			for (int i = 0; i < graph.size(); i++) {
				pregraph[i] += cos((M_PI * i * .5 * it / (factorDepth)));
			}
			if (checkCloser(graph, pregraph, sequence, factorDepth)) {
				for (int i = 0; i < graph.size(); i++) {
					graph[i] += cos((M_PI * i * .5 * it / (factorDepth)));
				}
				frequencies[ik].push_back(it / (factorDepth));
			}
			if (!checkCloser(graph, pregraph, sequence, factorDepth)) {
				for (int i = 0; i < graph.size(); i++) {
					pregraph[i] -= cos((M_PI * i * .5 * it / (factorDepth)));
				}
			}
		}
		for (int it = 0; it < factorDepth; it++) {
			for (int i = 0; i < graph.size(); i++) {
				pregraph[i] += cos((M_PI * .5 * i * it / (factorDepth)));
			}
			if (checkCloser(graph, pregraph, sequence, factorDepth)) {
				for (int i = 0; i < graph.size(); i++) {
					graph[i] += cos((M_PI * i * .5 * it / (factorDepth)));
				}
				frequencies[ik].push_back(it / (factorDepth));
			}
			if (!checkCloser(graph, pregraph, sequence, factorDepth)) {
				for (int i = 0; i < graph.size(); i++) {
					pregraph[i] -= cos((M_PI * i * .5 * it / (factorDepth)));
				}
			}
		}
		for (int it = 0; it < factorDepth; it++) {
			for (int i = 0; i < graph.size(); i++) {
				pregraph[i] += cos((M_PI * i * .5 * it / (factorDepth)));
			}
			if (checkCloser(graph, pregraph, sequence, factorDepth)) {
				for (int i = 0; i < graph.size(); i++) {
					graph[i] += cos((M_PI * i * .5 * it / (factorDepth)));
				}
				frequencies[ik].push_back(it / (factorDepth));
			}
			if (!checkCloser(graph, pregraph, sequence, factorDepth)) {
				for (int i = 0; i < graph.size(); i++) {
					pregraph[i] -= cos((M_PI * .5 * i * it / (factorDepth)));
				}
			}
		}
		for (int it = 0; it < factorDepth; it++) {
			for (int i = 0; i < graph.size(); i++) {
				pregraph[i] += cos((M_PI * i * .5 * it / (factorDepth)));
			}
			if (checkCloser(graph, pregraph, sequence, factorDepth)) {
				for (int i = 0; i < graph.size(); i++) {
					graph[i] += cos((M_PI * i * .5 * it / (factorDepth)));
				}
				frequencies[ik].push_back(it / (factorDepth));
			}
			if (!checkCloser(graph, pregraph, sequence, factorDepth)) {
				for (int i = 0; i < graph.size(); i++) {
					pregraph[i] -= cos((M_PI * i * .5 * it / (factorDepth)));
				}
			}
		}
	}

	std::cout << "\n \n \n Done. \n \n Roots sequence cosine frequencies: \n";
	for (int i = 0; i < frequencies[0].size(); i++) {
		std::cout << frequencies[0][i] << ", ";
	}
	std::cout << "\n Quotient sequence cosine frequencies: \n";
	for (int i = 0; i < frequencies[0].size(); i++) {
		std::cout << frequencies[1][i] << ", ";
	}
	std::cout << "\n Difference sequenc cosine frequencies: \n";
	for (int i = 0; i < frequencies[0].size(); i++) {
		std::cout << frequencies[2][i] << ", ";
	}
	std::cout << "\n Addend sequence cosine frequencies: \n";
	for (int i = 0; i < frequencies[0].size(); i++) {
		std::cout << frequencies[3][i] << ", ";
	}
	std::cout << "\n Product sequence cosine frequencies: \n";
	for (int i = 0; i < frequencies[0].size(); i++) {
		std::cout << frequencies[4][i] << ", ";
	}
	std::cout << "\n Powers sequence cosine frequencies: \n";
	for (int i = 0; i < frequencies[0].size(); i++) {
		std::cout << frequencies[5][i] << ", ";
	}

	std::cout << "\n \n \n Done sequencing frequency approximation scheduling. Sequence of frequencies (depth " << factorDepth << "). \n";
	
}



bool powcheckCloser(std::vector<double>graph, std::vector<double>pregraph, std::vector<double>sequence, double factorDepth) {
	std::vector<std::vector<bool>>nearingAnswerVector(2);
	nearingAnswerVector.resize(0);
	nearingAnswerVector.resize(2);
	nearingAnswerVector[0].resize(0);
	nearingAnswerVector[1].resize(0);
	for (int im = 0; im < sequence.size(); im++) {
		if (abs(graph.at((sequence[im] - 1) * factorDepth) - (sequence[im]-1))   <= abs(pregraph.at((sequence[im] - 1) * factorDepth - (sequence[im] - 1)))) {
			for (int i = 0; i < factorDepth; i++) {
				nearingAnswerVector[0].push_back(true);
			}
		}
		if (abs(graph.at((sequence[im] - 1) * factorDepth) - (sequence[im] - 1)) > abs(pregraph.at((sequence[im] - 1) * factorDepth) - (sequence[im] - 1))) {
			for (int i = 0; i < factorDepth; i++) {
				nearingAnswerVector[1].push_back(true);
			}
		}
	}

	for (int iz = 0; iz <graph.size() && iz<pregraph.size(); iz++) {
		if (round(iz / factorDepth) < sequence.size()) {
			if (abs(pregraph[iz] - sequence[round(iz / factorDepth)]) < abs(graph[iz] - sequence[round(iz / factorDepth)])) {
				nearingAnswerVector[0].push_back(true);
			}
			if (abs(graph[iz] - sequence[round(iz / factorDepth)]) <= abs(pregraph[iz] - sequence[round(iz / factorDepth)])) {
				nearingAnswerVector[1].push_back(true);
			}
			for (int iq = 0; iq < factorDepth; iq++) {
				if (iz*factorDepth+iq == sequence[round(iz / factorDepth)]) {
					if (abs(pregraph[iz] - sequence[round(iz / factorDepth)]) > abs(graph[iz] - sequence[round(iz / factorDepth)])) {
						if (nearingAnswerVector[0].size() > 0) {
							nearingAnswerVector[0].resize(nearingAnswerVector[0].size() - 1);
						}
						
					}
					if (abs(graph[iz] - sequence[round(iz / factorDepth)]) <= abs(pregraph[iz] - sequence[round(iz / factorDepth)])) {
						if (nearingAnswerVector[1].size() > 0) {
							nearingAnswerVector[1].resize(nearingAnswerVector[1].size() - 1);
						}
					}
				}
			}
		}
		
	}
	if (nearingAnswerVector[0].size() > nearingAnswerVector[1].size()) {
		//std::cout << "\n found a factor...";
		return true;
		
	}
	if (nearingAnswerVector[0].size() <= nearingAnswerVector[1].size()) {
		//std::cout << "\n did not find a factor...";
		return false;
	}

}


void curveSeq(std::vector<std::vector<std::vector<double>>>seq, double factorDepth) {
	double depth;
	std::vector<double>sequence;
	std::vector<double>graph;
	std::vector<double>pregraph;
	std::vector<std::vector<double>>factors;


	factors.resize(0);
	factors.resize(6);

	std::cout << "\n Searching for combinations of curve power factors to fit sequence, scanning range six times for patterning.";
	for (int ie = 0; ie < 6; ie++) {
		std::cout << "\n " << ie << "Set processing powers, \n";
		sequence.resize(0);
		depth = seq[ie][seq[ie].size() - 1].size();
		for (int it = 0; it < seq[ie].size(); it++) {
			sequence.push_back(seq[ie][it].size());
		}
		graph.resize(0);
		pregraph.resize(0);
		graph.resize(factorDepth * (depth + 1));
		pregraph.resize(factorDepth * (depth + 1));
		for (int it = 0; it < factorDepth; it++) {
			for (int i = 0; i < graph.size(); i++) {
				pregraph[i] += pow(i / factorDepth, it / (factorDepth / 2.0));
			}
			if (powcheckCloser(graph, pregraph, sequence, factorDepth)) {
				for (int i = 0; i < graph.size(); i++) {
					graph[i] += pow(i / factorDepth, it / (factorDepth / 2.0));
				}
				factors[ie].push_back(it / (factorDepth / 2.0));
			}
			if (!powcheckCloser(graph, pregraph, sequence, factorDepth)) {
				for (int i = 0; i < graph.size(); i++) {
					pregraph[i] -= pow(i / factorDepth, it / (factorDepth / 2.0));
				}
			}
		}
		for (int it = 0; it < factorDepth; it++) {
			for (int i = 0; i < graph.size(); i++) {
				pregraph[i] += pow(i / factorDepth, it / (factorDepth / 2.0));
			}
			if (powcheckCloser(graph, pregraph, sequence, factorDepth)) {
				for (int i = 0; i < graph.size(); i++) {
					graph[i] += pow(i / factorDepth, it / (factorDepth / 2.0));
				}
				factors[ie].push_back(it / (factorDepth / 2.0));
			}
			if (!powcheckCloser(graph, pregraph, sequence, factorDepth)) {
				for (int i = 0; i < graph.size(); i++) {
					pregraph[i] -= pow(i / factorDepth, it / (factorDepth / 2.0));
				}
			}
		}
		for (int it = 0; it < factorDepth; it++) {
			for (int i = 0; i < graph.size(); i++) {
				pregraph[i] += pow(i / factorDepth, it / (factorDepth / 2.0));
			}
			if (powcheckCloser(graph, pregraph, sequence, factorDepth)) {
				for (int i = 0; i < graph.size(); i++) {
					graph[i] += pow(i / factorDepth, it / (factorDepth / 2.0));
				}
				factors[ie].push_back(it / (factorDepth / 2.0));
			}
			if (!powcheckCloser(graph, pregraph, sequence, factorDepth)) {
				for (int i = 0; i < graph.size(); i++) {
					pregraph[i] -= pow(i / factorDepth, it / (factorDepth / 2.0));
				}
			}
		}
		for (int it = 0; it < factorDepth; it++) {
			for (int i = 0; i < graph.size(); i++) {
				pregraph[i] += pow(i / factorDepth, it / (factorDepth / 2.0));
			}
			if (powcheckCloser(graph, pregraph, sequence, factorDepth)) {
				for (int i = 0; i < graph.size(); i++) {
					graph[i] += pow(i / factorDepth, it / (factorDepth / 2.0));
				}
				factors[ie].push_back(it / (factorDepth / 2.0));
			}
			if (!powcheckCloser(graph, pregraph, sequence, factorDepth)) {
				for (int i = 0; i < graph.size(); i++) {
					pregraph[i] -= pow(i / factorDepth, it / (factorDepth / 2.0));
				}
			}
		}
		for (int it = 0; it < factorDepth; it++) {
			for (int i = 0; i < graph.size(); i++) {
				pregraph[i] += pow(i / factorDepth, it / (factorDepth / 2.0));
			}
			if (powcheckCloser(graph, pregraph, sequence, factorDepth)) {
				for (int i = 0; i < graph.size(); i++) {
					graph[i] += pow(i / factorDepth, it / (factorDepth / 2.0));
				}
				factors[ie].push_back(it / (factorDepth / 2.0));
			}
			if (!powcheckCloser(graph, pregraph, sequence, factorDepth)) {
				for (int i = 0; i < graph.size(); i++) {
					pregraph[i] -= pow(i / factorDepth, it / (factorDepth / 2.0));
				}
			}
		}
		for (int it = 0; it < factorDepth; it++) {
			for (int i = 0; i < graph.size(); i++) {
				pregraph[i] += pow(i / factorDepth, it / (factorDepth / 2.0));
			}
			if (powcheckCloser(graph, pregraph, sequence, factorDepth)) {
				for (int i = 0; i < graph.size(); i++) {
					graph[i] += pow(i / factorDepth, it / (factorDepth / 2.0));
				}
				factors[ie].push_back(it / (factorDepth / 2.0));
			}
			if (!powcheckCloser(graph, pregraph, sequence, factorDepth)) {
				for (int i = 0; i < graph.size(); i++) {
					pregraph[i] -= pow(i / factorDepth, it / (factorDepth / 2.0));
				}
			}
		}
		
	}


	std::cout << "\n \n \n Done. \n \n Roots sequence curve factors: \n";
	for (int i = 0; i < factors[0].size(); i++) {
		std::cout << factors[0][i] << ", ";
	}
	std::cout << "\n Quotient sequence curve factors: \n";
	for (int i = 0; i < factors[0].size(); i++) {
		std::cout << factors[1][i] << ", ";
	}
	std::cout << "\n Difference sequence curve factors: \n";
	for (int i = 0; i < factors[0].size(); i++) {
		std::cout << factors[2][i] << ", ";
	}
	std::cout << "\n Addend sequence curve factors: \n";
	for (int i = 0; i < factors[0].size(); i++) {
		std::cout << factors[3][i] << ", ";
	}
	std::cout << "\n Product sequence curve factors: \n";
	for (int i = 0; i < factors[0].size(); i++) {
		std::cout << factors[4][i] << ", ";
	}
	std::cout << "\n Powers sequence curve factors: \n";
	for (int i = 0; i < factors[0].size(); i++) {
		std::cout << factors[5][i] << ", ";
	}

	std::cout << "\n \n \n Done sequencing powers approximation scheduling. Sequence of factors (depth " << factorDepth << "). \n";


}





int main()
{
	double depth;
	double factorizationDepth;
	bool done=false;


	while (!done) {

		std::cout << "     Welcome to UniqueSolution, a problem solver: \n" <<
			" The Number of Unique Identities Created by the Six Major Operations on Any Two Integers 0-X \n" <<
			" Enter Depth of the Six Functions sequences to be indexed: ";
		std::cin >> depth;
		std::cout << "\n Processing depth: " << depth << ".";
		std::cout << "\n And, for after the sequences are loaded:"<<
			"\n Factorization hunt depth (number of factors to consider before rerouting (suggest minimum " << depth << ")) : ";
		std::cin >> factorizationDepth;
		calcSeq(depth);
		

		std::cout << "\n \n OK......  Cosine Frequencies Solutions Solver run. \n";
		modulaSeq(uniqueValues, factorizationDepth);

		std::cout << "\n \n OK......  Relations Curve Solutions Solver run. \n";
		curveSeq(uniqueValues, factorizationDepth);

		std::cout << "\n done?(0/1): ";
		std::cin >> done;
		std::cout << "\n";
	}
	std::cin >> done;
	return !done;
}


//root,divide,subtract,add,multiply,pow
