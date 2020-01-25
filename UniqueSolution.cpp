
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
								if (uniqueValues[i][it][io] == 1.0) {
									found=true;
								}
							}
						}
						if (!found) {
							
							if (!ii == 0) {
								uniqueValues[i][it].push_back(std::pow(iu, 1.0 / ii));
							}
							if (ii == 0) {
								uniqueValues[i][it].push_back(1.0);
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
								if (uniqueValues[i][it][io] == (depth*2)) {
									found=true;
								}
							}
						}
						if (!found) {

							if (!ii == 0) {
								uniqueValues[i][it].push_back(iu/ii);
							}
							if (ii == 0) {
								uniqueValues[i][it].push_back((depth*2));
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
	for (int it = 0; it < sequence.size(); it++) {
		if (abs(graph.at(sequence[it]*factorDepth))<= abs(pregraph.at(sequence[it] * factorDepth))) {
			for (int i = 0; i < factorDepth; i++) {
				nearingAnswerVector[0].push_back(true);
			}
		}
		if (abs(graph.at(sequence[it] * factorDepth)) > abs(pregraph.at(sequence[it] * factorDepth))) {
			for (int i = 0; i < factorDepth; i++) {
				nearingAnswerVector[1].push_back(true);
			}
		}
	}
	
	for (int i = 0; i < (sequence.size()-1) * factorDepth; i++) {
		if (abs(pregraph[i]) < abs(graph[i])) {
			nearingAnswerVector[0].push_back(true);
		}
		if (abs(graph[i]) <= abs(pregraph[i])) {
			nearingAnswerVector[1].push_back(true);
		}
		for (int it = 0; it < sequence.size();it++) {
			if (i / factorDepth == sequence[it]) {
				if (abs(pregraph[i]) > abs(graph[i])) {
					nearingAnswerVector[0].resize(nearingAnswerVector[0].size()-1);
				}
				if (abs(graph[i]) <= abs(pregraph[i])) {
					nearingAnswerVector[1].resize(nearingAnswerVector[1].size() - 1);
				}
			}
		}
	}
	if (nearingAnswerVector[0].size() > nearingAnswerVector[1].size()) {
		return true;
	}
	if (nearingAnswerVector[0].size() <= nearingAnswerVector[1].size()) {
		return false;
	}
	
}



void modulaSeq(std::vector<std::vector<double>>seq, double depth, double factorDepth) {

	std::vector<double>sequence;
	std::vector<double>graph;
	std::vector<double>pregraph;
	std::vector<std::vector<double>>frequencies;
	
	
	frequencies.resize(0);
	frequencies.resize(6);
	graph.resize(factorDepth * depth + 1);
	pregraph.resize(factorDepth * depth + 1);

	for (int i = 0; i < seq.size(); i++) {
		sequence.push_back(seq[i].size());
	}
	std::cout << "\n Searching for combinations of cosine modulations to fit sequence.";
	for (int i = 0; i < graph.size(); i++) {
		graph.resize(0);
		pregraph.resize(0);
		graph.resize(factorDepth * depth + 1);
		pregraph.resize(factorDepth * depth + 1);
		for (int it = 0; it < factorDepth; it++) {
			pregraph[i] += cos((M_PI * i * it / factorDepth) / sequence[sequence.size() - 1.0]);
			if (checkCloser(graph, pregraph, sequence, factorDepth)) {
				graph[i] += cos((M_PI * i * it / factorDepth) / sequence[sequence.size() - 1.0]);
				frequencies[i].push_back(it);
			}
			if (!checkCloser(graph, pregraph, sequence, factorDepth)) {
				pregraph[i] -= cos((M_PI * i * it / factorDepth) / sequence[sequence.size() - 1.0]);
			}
			for (int iu = 0; iu < factorDepth; iu++) {
				pregraph[i] += cos((M_PI * i * iu / factorDepth) / sequence[sequence.size() - 1.0]);
				if (checkCloser(graph, pregraph, sequence, factorDepth)) {
					graph[i] += cos((M_PI * i * iu / factorDepth) / sequence[sequence.size() - 1.0]);
					frequencies[i].push_back(iu);
				}
				if (!checkCloser(graph, pregraph, sequence, factorDepth)) {
					pregraph[i] -= cos((M_PI * i * iu / factorDepth) / sequence[sequence.size() - 1.0]);
				}
				for (int ii = 0; ii < factorDepth; ii++) {
					pregraph[i] += cos((M_PI * i * ii / factorDepth) / sequence[sequence.size() - 1.0]);
					if (checkCloser(graph, pregraph, sequence, factorDepth)) {
						graph[i] += cos((M_PI * i * ii / factorDepth) / sequence[sequence.size() - 1.0]);
						frequencies[i].push_back(ii);
					}
					if (!checkCloser(graph, pregraph, sequence, factorDepth)) {
						pregraph[i] -= cos((M_PI * i * ii / factorDepth) / sequence[sequence.size() - 1.0]);
					}
					for (int io = 0; io < factorDepth; io++) {
						pregraph[i] += cos((M_PI * i * io / factorDepth) / sequence[sequence.size() - 1.0]);
						if (checkCloser(graph, pregraph, sequence, factorDepth)) {
							graph[i] += cos((M_PI * i * io / factorDepth) / sequence[sequence.size() - 1.0]);
							frequencies[i].push_back(io);
						}
						if (!checkCloser(graph, pregraph, sequence, factorDepth)) {
							pregraph[i] -= cos((M_PI * i * io / factorDepth) / sequence[sequence.size() - 1.0]);
						}
						for (int ip = 0; ip < factorDepth; ip++) {
							pregraph[i] += cos((M_PI * i * ip / factorDepth) / sequence[sequence.size() - 1.0]);
							if (checkCloser(graph, pregraph, sequence, factorDepth)) {
								graph[i] += cos((M_PI * i * ip / factorDepth) / sequence[sequence.size() - 1.0]);
								frequencies[i].push_back(ip);
							}
							if (!checkCloser(graph, pregraph, sequence, factorDepth)) {
								pregraph[i] -= cos((M_PI * i * ip / factorDepth) / sequence[sequence.size() - 1.0]);
							}
							for (int ig = 0; ig < factorDepth; ig++) {
								pregraph[i] += cos((M_PI * i * ig / factorDepth) / sequence[sequence.size() - 1.0]);
								if (checkCloser(graph, pregraph, sequence, factorDepth)) {
									graph[i] += cos((M_PI * i * ig / factorDepth) / sequence[sequence.size() - 1.0]);
									frequencies[i].push_back(ig);
								}
								if (!checkCloser(graph, pregraph, sequence, factorDepth)) {
									pregraph[i] -= cos((M_PI * i * ig / factorDepth) / sequence[sequence.size() - 1.0]);
								}
							}
						}
					}
				}
			}
			std::cout << ".";
		}
		std::cout << "\n Next set processing cosines, \n";
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

bool seqcheckCloser(std::vector<double>graph, std::vector<double>pregraph, std::vector<double>sequence, double factorDepth) {
	std::vector<std::vector<bool>>nearingAnswerVector(2);
	nearingAnswerVector.resize(0);
	nearingAnswerVector.resize(2);
	for (int it = 0; it < sequence.size(); it++) {
		if (abs(graph.at(sequence[it] * factorDepth) - sequence[it])   <= abs(pregraph.at(sequence[it] * factorDepth - sequence[it]))) {
			for (int i = 0; i < factorDepth; i++) {
				nearingAnswerVector[0].push_back(true);
			}
		}
		if (abs(graph.at(sequence[it] * factorDepth) - sequence[it]) > abs(pregraph.at(sequence[it] * factorDepth) - sequence[it])) {
			for (int i = 0; i < factorDepth; i++) {
				nearingAnswerVector[1].push_back(true);
			}
		}
	}

	for (int i = 0; i < (sequence.size() - 1) * factorDepth; i++) {
		if (abs(pregraph[i] - sequence[round(i/factorDepth)]) < abs(graph[i] - sequence[round(i / factorDepth)])) {
			nearingAnswerVector[0].push_back(true);
		}
		if (abs(graph[i] - sequence[round(i / factorDepth)]) <= abs(pregraph[i] - sequence[round(i / factorDepth)])) {
			nearingAnswerVector[1].push_back(true);
		}
		for (int it = 0; it < sequence.size(); it++) {
			if (i / factorDepth == sequence[it]) {
				if (abs(pregraph[i] - sequence[round(i / factorDepth)]) > abs(graph[i] - sequence[round(i / factorDepth)])) {
					nearingAnswerVector[0].resize(nearingAnswerVector[0].size() - 1);
				}
				if (abs(graph[i] - sequence[round(i / factorDepth)]) <= abs(pregraph[i] - sequence[round(i / factorDepth)])) {
					nearingAnswerVector[1].resize(nearingAnswerVector[1].size() - 1);
				}
			}
		}
	}
	if (nearingAnswerVector[0].size() > nearingAnswerVector[1].size()) {
		return true;
	}
	if (nearingAnswerVector[0].size() <= nearingAnswerVector[1].size()) {
		return false;
	}

}
void curveSeq(std::vector<std::vector<double>>seq, double depth, double factorDepth) {

	std::vector<double>sequence;
	std::vector<double>graph;
	std::vector<double>pregraph;
	std::vector<std::vector<double>>factors;


	factors.resize(0);
	factors.resize(6);
	graph.resize(factorDepth * depth + 1);
	pregraph.resize(factorDepth * depth + 1);

	for (int i = 0; i < seq.size(); i++) {
		sequence.push_back(seq[i].size());
	}
	std::cout << "\n Searching for combinations of curve power factors to fit sequence.";
	for (int i = 0; i < graph.size(); i++) {
		graph.resize(0);
		pregraph.resize(0);
		graph.resize(factorDepth * depth + 1);
		pregraph.resize(factorDepth * depth + 1);
		for (int it = 0; it < factorDepth; it++) {
			pregraph[i] += pow(i/factorDepth,it /(factorDepth/2.0));
			if (checkCloser(graph, pregraph, sequence, factorDepth)) {
				graph[i] += pow(i / factorDepth, it / (factorDepth / 2.0));
				factors[i].push_back(it);
			}
			if (!checkCloser(graph, pregraph, sequence, factorDepth)) {
				pregraph[i] -= pow(i / factorDepth, it / (factorDepth / 2.0));
			}
			for (int iu = 0; iu < factorDepth; iu++) {
				pregraph[i] += pow(i / factorDepth, iu / (factorDepth / 2.0));
				if (checkCloser(graph, pregraph, sequence, factorDepth)) {
					graph[i] += pow(i / factorDepth, iu / (factorDepth / 2.0));
					factors[i].push_back(iu);
				}
				if (!checkCloser(graph, pregraph, sequence, factorDepth)) {
					pregraph[i] -= pow(i / factorDepth, iu / (factorDepth / 2.0));
				}
				for (int ii = 0; ii < factorDepth; ii++) {
					pregraph[i] += pow(i / factorDepth, ii / (factorDepth / 2.0));
					if (checkCloser(graph, pregraph, sequence, factorDepth)) {
						graph[i] += pow(i / factorDepth, ii / (factorDepth / 2.0));
						factors[i].push_back(ii);
					}
					if (!checkCloser(graph, pregraph, sequence, factorDepth)) {
						pregraph[i] -= pow(i / factorDepth, ii / (factorDepth / 2.0));
					}
					for (int io = 0; io < factorDepth; io++) {
						pregraph[i] += pow(i / factorDepth, io / (factorDepth / 2.0));
						if (checkCloser(graph, pregraph, sequence, factorDepth)) {
							graph[i] += pow(i / factorDepth, io / (factorDepth / 2.0));
							factors[i].push_back(io);
						}
						if (!checkCloser(graph, pregraph, sequence, factorDepth)) {
							pregraph[i] -= pow(i / factorDepth, io / (factorDepth / 2.0));
						}
						for (int ip = 0; ip < factorDepth; ip++) {
							pregraph[i] += pow(i / factorDepth, ip / (factorDepth / 2.0));
							if (checkCloser(graph, pregraph, sequence, factorDepth)) {
								graph[i] += pow(i / factorDepth, ip / (factorDepth / 2.0));
								factors[i].push_back(ip);
							}
							if (!checkCloser(graph, pregraph, sequence, factorDepth)) {
								pregraph[i] -= pow(i / factorDepth, ip / (factorDepth / 2.0));
							}
							for (int ig = 0; ig < factorDepth; ig++) {
								pregraph[i] += pow(i / factorDepth, ig / (factorDepth / 2.0));
								if (checkCloser(graph, pregraph, sequence, factorDepth)) {
									graph[i] += pow(i / factorDepth, ig / (factorDepth / 2.0));
									factors[i].push_back(ig);
								}
								if (!checkCloser(graph, pregraph, sequence, factorDepth)) {
									pregraph[i] -= pow(i / factorDepth, ig / (factorDepth / 2.0));
								}
							}
						}
					}
				}
			}
			std::cout << ".";
		}
		std::cout << "\n Next set processing powers, \n";
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

		std::cout << "Welcome to UniqueSolution, a problem solver: \n" <<
			"The Number of Unique Identities Created by the Six Major Operations on Any Two Integers 0-X \n" <<
			"Enter Depth of the Six Functions sequences to be indexed: ";
		std::cin >> depth;
		std::cout << "\n Processing depth: " << depth << ".";
		calcSeq(depth);
		std::cout << "\n \n \n Factorization hunt depth (number of factors to consider before rerouting (suggest minimum 32)) : ";
		std::cin >> factorizationDepth;

		std::cout << "\n \n OK......  Cosine Frequencies Solutions Solver run. \n";
		for (int i = 0; i < 6; i++) {
			modulaSeq(uniqueValues[i], uniqueValues[i][uniqueValues[i].size()-1.0].size(), factorizationDepth);
		}

		std::cout << "\n \n OK......  Relations Curve Solutions Solver run. \n";
		for (int i = 0; i < 6; i++) {
			curveSeq(uniqueValues[i], uniqueValues[i][uniqueValues[i].size()-1.0].size(), factorizationDepth);
		}

		std::cout << "\n done?(0/1): ";
		std::cin >> done;
		std::cout << "\n";
	}
	std::cin >> done;
	return !done;
}


//root,divide,subtract,add,multiply,pow
