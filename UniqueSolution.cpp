/// UniqueSolution.cpp : The Number of Unique Identities Created by the Six Major Operations on Any Two Integers 0-X,
//     Function Finder
//
#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <cmath>
#include<math.h>
#include <vector>


std::vector<std::vector<std::vector<double>>> uniqueValues(6);
int cycles;

void calcSeq(double depth) {
	for (double it = 0; it < 6; it++) {
		uniqueValues[it].resize(depth);
	}
	for (double i = 0; i < 6; i++) {
		for (double it = 0; it < depth; it++) {
			for (double iu = 0; iu <= it; iu++) {
				for (double ii = 0; ii <= it; ii++) {
					if (i == 0) {
						bool found = false;
						for (int io = 0; io < uniqueValues[i][it].size(); io++) {
							if (!ii == 0) {
								if (uniqueValues[i][it][io] == std::pow(iu, 1.0 / ii)) {
									found = true;
								}
							}
							if (ii == 0) {
								if (iu != 0 && uniqueValues[i][it][io] == (M_PI)) {
									found = true;
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
									found = true;
								}
							}
							if (ii == 0) {
								if (iu <= 0) {
									if (uniqueValues[i][it][io] == (0)) {
										found = true;
									}
								}
								else if (uniqueValues[i][it][io] == (DBL_MAX)) {
									found = true;
								}
							}
						}
						if (!found) {

							if (!ii == 0) {
								uniqueValues[i][it].push_back(iu / ii);
							}

							if (ii == 0) {
								if (iu != 0) {
									uniqueValues[i][it].push_back((DBL_MAX));
								}
								if (iu == 0) {
									uniqueValues[i][it].push_back(0);
								}
							}
						}
					}
					if (i == 2) {
						bool found = false;
						for (int io = 0; io < uniqueValues[i][it].size(); io++) {
							if (uniqueValues[i][it][io] == (iu - ii)) {
								found = true;
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
								found = true;
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
								found = true;
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
								found = true;
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
	for (int i = 0; i < depth; i++) {
		std::cout << uniqueValues[0][i].size() << ", ";
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
std::cout << uniqueValues[5][i].size() << ", ";
	}//root,divide,subtract,add,multiply,pow
}//root,divide,subtract,add,multiply,pow


bool findInList(std::vector<double>list, double num) {
	for (int m = 0; m < list.size(); m++) {
		if (abs(list[m]) == abs(num)) {
			return true;
		}

	}
	return false;
}



int checkCloser(std::vector<double>pregrapha, std::vector<double>pregraphs, std::vector<double>invpregrapha, 
	std::vector<double>invpregraphs, std::vector<double>graph, std::vector<double>sequence, double factorDepth) {

	std::vector<std::vector<bool>>nearingAnswerVector(2);
	nearingAnswerVector.resize(0);
	nearingAnswerVector.resize(8);
	nearingAnswerVector[0].resize(0);
	nearingAnswerVector[1].resize(0);
	nearingAnswerVector[2].resize(0);
	nearingAnswerVector[3].resize(0);
	nearingAnswerVector[4].resize(0);
	nearingAnswerVector[5].resize(0);
	nearingAnswerVector[6].resize(0);
	nearingAnswerVector[7].resize(0);
	for (int ic = 0; ic < sequence.size(); ic++) {
		if (abs(graph.at(sequence[ic] * factorDepth)) >= abs(pregrapha.at(sequence[ic] * factorDepth))) {
			if (pregrapha.at(sequence[ic] * factorDepth) == 0) {
				for (int f = 0; f < (factorDepth*2) * sequence[sequence.size()-1]/sequence.size(); f++) {
					nearingAnswerVector[0].push_back(true);
				}
			}
			else {
				for (int f = 0; f < (((factorDepth)* sequence[sequence.size() - 1]) / (sequence.size()))*abs(abs(graph[ic*factorDepth])- abs((pregrapha[ic * factorDepth])));f++) {
					nearingAnswerVector[0].push_back(true);
				}
			}
			
		}
		if (abs(graph.at(sequence[ic] * factorDepth)) < abs(pregrapha.at(sequence[ic] * factorDepth))) {
			for (int f = 0; f < (((factorDepth*2)* sequence[sequence.size() - 1]) / (sequence.size())) * abs((pregrapha[ic * factorDepth])); f++) {
				nearingAnswerVector[1].push_back(true);
			}
		}


		if (abs(graph.at(sequence[ic] * factorDepth)) >= abs(pregraphs.at(sequence[ic] * factorDepth))) {
			if (pregraphs.at(sequence[ic] * factorDepth) == 0) {
				for (int f = 0; f < (factorDepth*2)* sequence[sequence.size() - 1] / sequence.size(); f++) {
					nearingAnswerVector[2].push_back(true);
				}
			}
			else {
				for (int f = 0; f < (((factorDepth*2)* sequence[sequence.size() - 1]) / (sequence.size())) * abs(abs((graph[ic * factorDepth])) - abs((pregraphs[ic * factorDepth]))); f++) {
					nearingAnswerVector[2].push_back(true);
				}
			}

		}
		if (abs(graph.at(sequence[ic] * factorDepth)) < abs(pregraphs.at(sequence[ic] * factorDepth))) {
			for (int f = 0; f < (((factorDepth *2)* sequence[sequence.size() - 1]) / (sequence.size())) * abs((pregraphs[ic * factorDepth])); f++) {
				nearingAnswerVector[3].push_back(true);
			}
		}



		if (abs(graph.at(sequence[ic] * factorDepth)) >= abs(invpregrapha.at(sequence[ic] * factorDepth))) {
			if (invpregrapha.at(sequence[ic] * factorDepth) == 0) {
				for (int f = 0; f < (factorDepth*2)* sequence[sequence.size() - 1] / sequence.size(); f++) {
					nearingAnswerVector[4].push_back(true);
				}
			}
			else {
				for (int f = 0; f < (((factorDepth)* sequence[sequence.size() - 1]) / (sequence.size())) * abs(abs((graph[ic * factorDepth])) - abs((invpregrapha[ic * factorDepth]))); f++) {
					nearingAnswerVector[4].push_back(true);
				}
			}

		}
		if (abs(graph.at(sequence[ic] * factorDepth)) < abs(invpregrapha.at(sequence[ic] * factorDepth))) {
			for (int f = 0; f < (((factorDepth*2 )* sequence[sequence.size() - 1]) / (sequence.size())) * abs((invpregrapha[ic * factorDepth])); f++) {
				nearingAnswerVector[5].push_back(true);
			}
		}


		if (abs(graph.at(sequence[ic] * factorDepth)) >= abs(invpregraphs.at(sequence[ic] * factorDepth))) {
			if (invpregraphs.at(sequence[ic] * factorDepth) == 0) {
				for (int f = 0; f < (factorDepth*2)* sequence[sequence.size() - 1] / sequence.size(); f++) {
					nearingAnswerVector[6].push_back(true);
				}
			}
			else {
				for (int f = 0; f < (((factorDepth)* sequence[sequence.size() - 1]) / (sequence.size())) * abs(abs((graph[ic * factorDepth])) - abs((invpregraphs[ic * factorDepth]))); f++) {
					nearingAnswerVector[6].push_back(true);
				}
			}

		}
		if (abs(graph.at(sequence[ic] * factorDepth)) < abs(invpregraphs.at(sequence[ic] * factorDepth))) {
			for (int f = 0; f < (((factorDepth *2)* sequence[sequence.size() - 1]) / (sequence.size())) * abs((invpregraphs[ic * factorDepth])); f++) {
				nearingAnswerVector[7].push_back(true);
			}

		}

	}

	for (int ib = 0; ib < graph.size(); ib++) {//maybe seqsize-1
		if (round(ib / factorDepth) != (ib / factorDepth)) {
			if (abs(pregrapha[ib]) >= abs(graph[ib])) {
				for(int f=0;f<abs(abs(pregrapha[ib])- abs(graph[ib]))*4;f++){
					nearingAnswerVector[0].push_back(true);
				}
			}
			if (abs(pregrapha[ib]) < abs(graph[ib])) {
				if (pregrapha[ib] != 0) {
					for (int f = 0; f < abs(abs((pregrapha[ib]) - abs(graph[ib])) * 4); f++) {
						nearingAnswerVector[1].push_back(true);
					}
				}
				else {
					for (int f = 0; f < abs(graph[ib] * 4); f++) {
						nearingAnswerVector[1].push_back(true);
					}
				}
			}
			if (abs(pregraphs[ib]) >= abs(graph[ib])) {
				for (int f = 0; f < abs(abs(pregraphs[ib]) - abs(graph[ib])) * 4; f++) {
					nearingAnswerVector[2].push_back(true);
				}
			}
			if (abs(pregraphs[ib]) < abs(graph[ib])) {
				if (pregraphs[ib] != 0) {
					for (int f = 0; f < abs(abs((pregraphs[ib])) - abs(graph[ib])) * 4; f++) {
						nearingAnswerVector[3].push_back(true);
					}
				}
				else {
					for (int f = 0; f < abs(graph[ib] * 4); f++) {
						nearingAnswerVector[3].push_back(true);
					}
				}
			}
			if (abs(invpregrapha[ib]) >= abs(graph[ib])) {
				for (int f = 0; f < abs(abs(invpregrapha[ib]) - abs(graph[ib])) * 4; f++) {
					nearingAnswerVector[4].push_back(true);
				}
			}
			if (abs(invpregrapha[ib]) < abs(graph[ib])) {
				if (invpregrapha[ib] != 0) {
					for (int f = 0; f < abs(abs((invpregrapha[ib])) - abs(graph[ib])) * 4; f++) {
						nearingAnswerVector[5].push_back(true);
					}
				}
				else {
					for (int f = 0; f < abs(graph[ib] * 4); f++) {
						nearingAnswerVector[5].push_back(true);
					}
				}
			}
			if (abs(invpregraphs[ib]) >= abs(graph[ib])) {
				for (int f = 0; f < abs(abs(invpregraphs[ib]) - abs(graph[ib])) * 4; f++) {
					nearingAnswerVector[6].push_back(true);
				}
			}
			if (abs(invpregraphs[ib]) < abs(graph[ib])) {
				if (invpregraphs[ib] != 0) {
					for (int f = 0; f < abs((abs(invpregraphs[ib])) - abs(graph[ib])) * 4; f++) {
						nearingAnswerVector[7].push_back(true);
					}
				}
				else {
					for (int f = 0; f < abs(graph[ib] * 4); f++) {
						nearingAnswerVector[7].push_back(true);
					}
				}
			}
		}
		
	}
	std::vector<double>getHigh;
	double highest;
	highest = 0;
	getHigh.resize(0);
	getHigh.resize(4);
	for (int f = 0; f < 4; f++) {
		getHigh[f] = nearingAnswerVector[2.0 * f].size() - nearingAnswerVector[(2.0 * f) + 1].size();
	}
	for (int f = 0; f < 4; f++) {
		if (getHigh[f] > highest) {
			highest = getHigh[f];
		}
	}
	for (int f = 0; f < 4; f++) {
		if (getHigh[f] == highest ) {
			if (highest > 1) {
				return f;
				break;
			}
			break;
		}
	}
	if (highest < 1) {
		return 4;
	}

}


//pregrapha,pregraphs,invpregrapha,invpregraphs
void modulaSeq(std::vector<std::vector<std::vector<double>>>seq, double factorDepth) {
	double depth;
	std::vector<double>sequence;
	std::vector<double>graph;
	std::vector<double>pregrapha;
	std::vector<double>pregraphs;
	std::vector<double>invpregrapha;
	std::vector<double>invpregraphs;
	std::vector<std::vector<double>>frequencies;

	frequencies.resize(0);
	frequencies.resize(6);
	std::cout << "\n Searching for combinations of cosine modulations to fit sequence, scanning range "<<cycles<<" times for patterning.";
	for (int ik = 0; ik < 6; ik++) {
		std::cout << "\n " << ik << "Set processing cosines, values found (will be listed later): \n";
		sequence.resize(0);
		depth = seq[ik][seq[ik].size() - 1].size();
		for (int it = 0; it < seq[ik].size(); it++) {
			sequence.push_back(seq[ik][it].size());
		}
		graph.resize(0);
		pregrapha.resize(0);
		pregraphs.resize(0);
		invpregrapha.resize(0);
		invpregraphs.resize(0);
		graph.resize(factorDepth * (depth +1));
		pregrapha.resize(factorDepth * (depth +1));
		pregraphs.resize(factorDepth * (depth +1));
		invpregrapha.resize(factorDepth * (depth + 1));
		invpregraphs.resize(factorDepth * (depth + 1));
		for (int z = 0; z < cycles; z++) {
			for (double it = 1; it < pow(factorDepth,2.0); it++) {//potentially undo that...
				//std::cout << "\nchecking case for desired value";
				for (int i = 0; i < graph.size(); i++) {
					pregrapha[i] += cos((M_PI * i * ((factorDepth / (it * factorDepth))) / factorDepth));
					invpregrapha[i] += cos((M_PI * i * ((it*factorDepth) / ((factorDepth))) / factorDepth));
					invpregraphs[i] -= cos((M_PI * i * ((it*factorDepth) / ((factorDepth))) / factorDepth));
					pregraphs[i] -= cos((M_PI * i * (((factorDepth / (it * factorDepth))) / factorDepth)));
				}
				switch (checkCloser(pregrapha, pregraphs, invpregrapha, invpregraphs, graph, sequence, factorDepth)) {
				case 0:
					//std::cout << "\nadding to graph";
					if (findInList(frequencies[ik], (factorDepth / (it * factorDepth))) == false) {
						std::cout << (factorDepth / (it * factorDepth))<<", ";
						for (int i = 0; i < graph.size(); i++) {
							pregraphs[i] += 2 * cos((M_PI * i * (factorDepth / (it * factorDepth))) / factorDepth);
							graph[i] += cos((M_PI * i * (factorDepth / (it * factorDepth))) / factorDepth);
							invpregrapha[i] -= cos((M_PI * i * ((it*factorDepth) / ((factorDepth))) / factorDepth));
							invpregraphs[i] += cos((M_PI * i * ((it*factorDepth) / ((factorDepth))) / factorDepth));
							invpregrapha[i] += cos((M_PI * i * (((factorDepth / (it * factorDepth))) / factorDepth)));
							invpregraphs[i] += cos((M_PI * i * (((factorDepth / (it * factorDepth))) / factorDepth)));
						}
						frequencies[ik].push_back((factorDepth / (it * factorDepth)));
					}
					else {
						for (int i = 0; i < graph.size(); i++) {
							pregraphs[i] += cos((M_PI * i * (factorDepth / (it * factorDepth))) / factorDepth);
							pregrapha[i] -= cos((M_PI * i * (factorDepth / (it * factorDepth))) / factorDepth);
							invpregrapha[i] -= cos((M_PI * i * ((it*factorDepth) / ((factorDepth))) / factorDepth));
							invpregraphs[i] += cos((M_PI * i * ((it*factorDepth) / ((factorDepth))) / factorDepth));
						}
					}
					break;
				case 1:
					//std::cout << "\nsubtracting to graph";
					if ( findInList(frequencies[ik], (factorDepth / (it * factorDepth))) == false) {
						std::cout << -(factorDepth / (it * factorDepth)) << ", ";
						for (int i = 0; i < graph.size(); i++) {
							pregrapha[i] -= 2 * cos((M_PI * i * (factorDepth / (it * factorDepth))) / factorDepth);
							graph[i] -= cos((M_PI * i * (factorDepth / (it * factorDepth))) / factorDepth);
							invpregrapha[i] -= cos((M_PI * i * ((it*factorDepth) / ((factorDepth))) / factorDepth));
							invpregraphs[i] += cos((M_PI * i * ((it*factorDepth) / ((factorDepth))) / factorDepth));
							invpregrapha[i] -= cos((M_PI * i * (((factorDepth / (it * factorDepth))) / factorDepth)));
							invpregraphs[i] -= cos((M_PI * i * (((factorDepth / (it * factorDepth))) / factorDepth)));
						}
						frequencies[ik].push_back(-(factorDepth / (it * factorDepth)));
					}
					else {
						for (int i = 0; i < graph.size() && findInList(frequencies[ik], (factorDepth / (it * factorDepth))) == false; i++) {
							pregraphs[i] += cos((M_PI * i * (factorDepth / (it * factorDepth))) / factorDepth);
							pregrapha[i] -= cos((M_PI * i * (factorDepth / (it * factorDepth))) / factorDepth);
							invpregrapha[i] -= cos((M_PI * i * ((it*factorDepth) / ((factorDepth))) / factorDepth));
							invpregraphs[i] += cos((M_PI * i * ((it*factorDepth) / ((factorDepth))) / factorDepth));
						}
					}
					break;
				case 2:
					//std::cout << "\nadding to graph";
					if (findInList(frequencies[ik], (it*factorDepth) /(factorDepth / (it * factorDepth))) == false) {
						std::cout << ((it * factorDepth) / ((factorDepth))) << ", ";
						for (int i = 0; i < graph.size(); i++) {
							pregrapha[i] -= cos((M_PI * i * (factorDepth / (it * factorDepth))) / factorDepth);
							pregraphs[i] += cos((M_PI * i * (((factorDepth / (it * factorDepth))) / factorDepth)));
							graph[i] += cos((M_PI * i * ((it*factorDepth) / ((factorDepth / (it * factorDepth))) / factorDepth)));
							invpregraphs[i] += 2 * cos((M_PI * i * ((it*factorDepth) / ((factorDepth))) / factorDepth));
							pregrapha[i] += cos((M_PI * i * ((it*factorDepth) / ((factorDepth))) / factorDepth));
							pregraphs[i] += cos((M_PI * i * ((it*factorDepth) / ((factorDepth))) / factorDepth));
						}
						frequencies[ik].push_back(((it* factorDepth) / ((factorDepth))));
					}
					else {
						for (int i = 0; i < graph.size(); i++) {
							pregraphs[i] += cos((M_PI * i * (factorDepth / (it * factorDepth))) / factorDepth);
							pregrapha[i] -= cos((M_PI * i * (factorDepth / (it * factorDepth))) / factorDepth);
							invpregrapha[i] -= cos((M_PI * i * ((it*factorDepth) / ((factorDepth))) / factorDepth));
							invpregraphs[i] += cos((M_PI * i * ((it*factorDepth) / ((factorDepth))) / factorDepth));
						}
					}
					break;
				case 3:
					//std::cout << "\nsubtracting to graph";
					if (findInList(frequencies[ik], (it*factorDepth) /(factorDepth / (it * factorDepth))) == false) {
						std::cout << -((it* factorDepth) / ((factorDepth))) << ", ";
						for (int i = 0; i < graph.size(); i++) {
							pregrapha[i] -= cos((M_PI * i * (factorDepth / (it * factorDepth))) / factorDepth);
							pregraphs[i] += cos((M_PI * i * (((factorDepth / (it * factorDepth))) / factorDepth)));
							graph[i] -= cos((M_PI * i * ((it*factorDepth) / ((factorDepth))) / factorDepth));
							invpregrapha[i] -= 2 * cos((M_PI * i * ((it*factorDepth) / ((factorDepth))) / factorDepth));
							pregrapha[i] -= cos((M_PI * i * ((it*factorDepth) / ((factorDepth))) / factorDepth));
							pregraphs[i] -= cos((M_PI * i * ((it*factorDepth) / ((factorDepth))) / factorDepth));
						}
						frequencies[ik].push_back(-((it * factorDepth) / ((factorDepth))));
					}
					else {
						for (int i = 0; i < graph.size(); i++) {
							pregraphs[i] += cos((M_PI * i * (factorDepth / (it * factorDepth))) / factorDepth);
							pregrapha[i] -= cos((M_PI * i * (factorDepth / (it * factorDepth))) / factorDepth);
							invpregrapha[i] -= cos((M_PI * i * ((it*factorDepth) / ((factorDepth))) / factorDepth));
							invpregraphs[i] += cos((M_PI * i * ((it*factorDepth)/ ((factorDepth / (it * factorDepth)))) / factorDepth));
						}
					}
					break;

				case 4:
					//std::cout << "\nno good value";
					if (graph[0] != pregrapha[0]) {
						for (int i = 0; i < graph.size(); i++) {
							pregraphs[i] += cos((M_PI * i * (factorDepth / (it * factorDepth))) / factorDepth);
							pregrapha[i] -= cos((M_PI * i * (factorDepth / (it * factorDepth))) / factorDepth);
							invpregrapha[i] -= cos((M_PI * i * ((it*factorDepth) / ((factorDepth))) / factorDepth));
							invpregraphs[i] += cos((M_PI * i * ((it*factorDepth) / ((factorDepth))) / factorDepth));
						}
					}
					break;
				}
			}
		}
	}

	std::cout << "\n \n \n Done. \n \n Roots sequence cosine frequencies: \n";
	for (int i = 0; i < frequencies[0].size(); i++) {
		std::cout << frequencies[0][i] << ", ";
	}
	std::cout << "\n Quotient sequence cosine frequencies: \n";
	for (int i = 0; i < frequencies[1].size(); i++) {
		std::cout << frequencies[1][i] << ", ";
	}
	std::cout << "\n Difference sequence cosine frequencies: \n";
	for (int i = 0; i < frequencies[2].size(); i++) {
		std::cout << frequencies[2][i] << ", ";
	}
	std::cout << "\n Addend sequence cosine frequencies: \n";
	for (int i = 0; i < frequencies[3].size(); i++) {
		std::cout << frequencies[3][i] << ", ";
	}
	std::cout << "\n Product sequence cosine frequencies: \n";
	for (int i = 0; i < frequencies[4].size(); i++) {
		std::cout << frequencies[4][i] << ", ";
	}
	std::cout << "\n Powers sequence cosine frequencies: \n";
	for (int i = 0; i < frequencies[5].size(); i++) {
		std::cout << frequencies[5][i] << ", ";
	}

	std::cout << "\n \n \n Done sequencing complex frequencies for sequence approximation. Sequence of frequencies (depth " << factorDepth << "). \n\n\n"
		"Now approximating singe-cosine frequecy modulation.";

}


void fmmodulaSeq(std::vector<std::vector<std::vector<double>>>seq, double factorDepth) {
	double depth;
	std::vector<double>sequence;
	std::vector<double>graph;
	std::vector<double>midgrapha;
	std::vector<double>pregrapha;
	std::vector<double>midgraphs;
	std::vector<double>pregraphs;
	std::vector<double>invmidgrapha;
	std::vector<double>invpregrapha;
	std::vector<double>invmidgraphs;
	std::vector<double>invpregraphs;
	std::vector<std::vector<double>>frequencies;
	frequencies.resize(0);
	frequencies.resize(12);
	std::cout << "\n Searching for combinations of cosine frequency modulations to fit sequence.";
	for (int ik = 0; ik < factorDepth*factorDepth; ik++) {
		for (int ik = 0; ik < 6; ik++) {
			std::cout << "\n " << ik << "Set processing cosines, values found (will be listed later): \n";
			sequence.resize(0);
			depth = seq[ik][seq[ik].size() - 1].size();
			for (int it = 0; it < seq[ik].size(); it++) {
				sequence.push_back(seq[ik][it].size());
			}
			graph.resize(0);
			pregraphs.resize(0);
			midgrapha.resize(0);
			pregrapha.resize(0);
			midgraphs.resize(0);
			invpregraphs.resize(0);
			invmidgrapha.resize(0);
			invpregrapha.resize(0);
			invmidgraphs.resize(0);
			graph.resize(factorDepth * (depth + 1));
			midgrapha.resize(factorDepth * (depth + 1));
			midgraphs.resize(factorDepth * (depth + 1));
			pregrapha.resize(factorDepth * (depth + 1));
			pregraphs.resize(factorDepth * (depth + 1));
			invmidgrapha.resize(factorDepth * (depth + 1));
			invmidgraphs.resize(factorDepth * (depth + 1));
			invpregrapha.resize(factorDepth * (depth + 1));
			invpregraphs.resize(factorDepth * (depth + 1));
			for (int z = 0; z < cycles; z++) {
				for (double it = 1; it < pow(factorDepth, 2.0); it++) {//potentially undo that...
					//std::cout << "\nchecking case for desired value";
					for (int i = 0; i < graph.size(); i++) {
						pregrapha[i] += cos(M_PI * i * (1 / it) / factorDepth);
						pregraphs[i] -= cos(M_PI * i * (1 / it) / factorDepth);
						invpregrapha[i] += cos(M_PI * i * (it / 1) / factorDepth);
						invpregraphs[i] -= cos(M_PI * i * (it / 1) / factorDepth);
					}
					for (int i = 0; i < graph.size(); i++) {
						midgrapha[i] = cos((M_PI * i * +pregrapha[i]) / factorDepth);
						midgraphs[i] = cos((M_PI * i * +pregraphs[i]) / factorDepth);
						invmidgrapha[i] = cos((M_PI * i * +invpregrapha[i]) / factorDepth);
						invmidgraphs[i] = cos((M_PI * i * +invpregraphs[i]) / factorDepth);

					}
					switch (checkCloser(midgrapha, midgraphs, invmidgrapha, invmidgraphs, graph, sequence, factorDepth)) {
					case 0:
						//std::cout << "\nadding to graph";
						if (findInList(frequencies[ik], (factorDepth / (it * factorDepth))) == false) {
							std::cout << (factorDepth / (it * factorDepth)) << ", ";
							for (int i = 0; i < graph.size(); i++) {
								pregraphs[i] += 2 * cos((M_PI * i * (factorDepth / (it * factorDepth))) / factorDepth);
								graph[i] += cos((M_PI * i * (factorDepth / (it * factorDepth))) / factorDepth);
								invpregrapha[i] -= cos((M_PI * i * ((it * factorDepth) / ((factorDepth))) / factorDepth));
								invpregraphs[i] += cos((M_PI * i * ((it * factorDepth) / ((factorDepth))) / factorDepth));
								invpregrapha[i] += cos((M_PI * i * (((factorDepth / (it * factorDepth))) / factorDepth)));
								invpregraphs[i] += cos((M_PI * i * (((factorDepth / (it * factorDepth))) / factorDepth)));
							}
							frequencies[ik].push_back((factorDepth / (it * factorDepth)));
						}
						else {
							for (int i = 0; i < graph.size(); i++) {
								pregraphs[i] += cos((M_PI * i * (factorDepth / (it * factorDepth))) / factorDepth);
								pregrapha[i] -= cos((M_PI * i * (factorDepth / (it * factorDepth))) / factorDepth);
								invpregrapha[i] -= cos((M_PI * i * ((it * factorDepth) / ((factorDepth))) / factorDepth));
								invpregraphs[i] += cos((M_PI * i * ((it * factorDepth) / ((factorDepth))) / factorDepth));
							}
						}
						break;
					case 1:
						//std::cout << "\nsubtracting to graph";
						if (findInList(frequencies[ik], (factorDepth / (it * factorDepth))) == false) {
							std::cout << -(factorDepth / (it * factorDepth)) << ", ";
							for (int i = 0; i < graph.size(); i++) {
								pregrapha[i] -= 2 * cos((M_PI * i * (factorDepth / (it * factorDepth))) / factorDepth);
								graph[i] -= cos((M_PI * i * (factorDepth / (it * factorDepth))) / factorDepth);
								invpregrapha[i] -= cos((M_PI * i * ((it * factorDepth) / ((factorDepth))) / factorDepth));
								invpregraphs[i] += cos((M_PI * i * ((it * factorDepth) / ((factorDepth))) / factorDepth));
								invpregrapha[i] -= cos((M_PI * i * (((factorDepth / (it * factorDepth))) / factorDepth)));
								invpregraphs[i] -= cos((M_PI * i * (((factorDepth / (it * factorDepth))) / factorDepth)));
							}
							frequencies[ik].push_back(-(factorDepth / (it * factorDepth)));
						}
						else {
							for (int i = 0; i < graph.size() && findInList(frequencies[ik], (factorDepth / (it * factorDepth))) == false; i++) {
								pregraphs[i] += cos((M_PI * i * (factorDepth / (it * factorDepth))) / factorDepth);
								pregrapha[i] -= cos((M_PI * i * (factorDepth / (it * factorDepth))) / factorDepth);
								invpregrapha[i] -= cos((M_PI * i * ((it * factorDepth) / ((factorDepth))) / factorDepth));
								invpregraphs[i] += cos((M_PI * i * ((it * factorDepth) / ((factorDepth))) / factorDepth));
							}
						}
						break;
					case 2:
						//std::cout << "\nadding to graph";
						if (findInList(frequencies[ik], (it * factorDepth) / (factorDepth / (it * factorDepth))) == false) {
							std::cout << ((it * factorDepth) / ((factorDepth))) << ", ";
							for (int i = 0; i < graph.size(); i++) {
								pregrapha[i] -= cos((M_PI * i * (factorDepth / (it * factorDepth))) / factorDepth);
								pregraphs[i] += cos((M_PI * i * (((factorDepth / (it * factorDepth))) / factorDepth)));
								graph[i] += cos((M_PI * i * ((it * factorDepth) / ((factorDepth / (it * factorDepth))) / factorDepth)));
								invpregraphs[i] += 2 * cos((M_PI * i * ((it * factorDepth) / ((factorDepth))) / factorDepth));
								pregrapha[i] += cos((M_PI * i * ((it * factorDepth) / ((factorDepth))) / factorDepth));
								pregraphs[i] += cos((M_PI * i * ((it * factorDepth) / ((factorDepth))) / factorDepth));
							}
							frequencies[ik].push_back(((it * factorDepth) / ((factorDepth))));
						}
						else {
							for (int i = 0; i < graph.size(); i++) {
								pregraphs[i] += cos((M_PI * i * (factorDepth / (it * factorDepth))) / factorDepth);
								pregrapha[i] -= cos((M_PI * i * (factorDepth / (it * factorDepth))) / factorDepth);
								invpregrapha[i] -= cos((M_PI * i * ((it * factorDepth) / ((factorDepth))) / factorDepth));
								invpregraphs[i] += cos((M_PI * i * ((it * factorDepth) / ((factorDepth))) / factorDepth));
							}
						}
						break;
					case 3:
						//std::cout << "\nsubtracting to graph";
						if (findInList(frequencies[ik], (it * factorDepth) / (factorDepth / (it * factorDepth))) == false) {
							std::cout << -((it * factorDepth) / ((factorDepth))) << ", ";
							for (int i = 0; i < graph.size(); i++) {
								pregrapha[i] -= cos((M_PI * i * (factorDepth / (it * factorDepth))) / factorDepth);
								pregraphs[i] += cos((M_PI * i * (((factorDepth / (it * factorDepth))) / factorDepth)));
								graph[i] -= cos((M_PI * i * ((it * factorDepth) / ((factorDepth))) / factorDepth));
								invpregrapha[i] -= 2 * cos((M_PI * i * ((it * factorDepth) / ((factorDepth))) / factorDepth));
								pregrapha[i] -= cos((M_PI * i * ((it * factorDepth) / ((factorDepth))) / factorDepth));
								pregraphs[i] -= cos((M_PI * i * ((it * factorDepth) / ((factorDepth))) / factorDepth));
							}
							frequencies[ik].push_back(-((it * factorDepth) / ((factorDepth))));
						}
						else {
							for (int i = 0; i < graph.size(); i++) {
								pregraphs[i] += cos((M_PI * i * (factorDepth / (it * factorDepth))) / factorDepth);
								pregrapha[i] -= cos((M_PI * i * (factorDepth / (it * factorDepth))) / factorDepth);
								invpregrapha[i] -= cos((M_PI * i * ((it * factorDepth) / ((factorDepth))) / factorDepth));
								invpregraphs[i] += cos((M_PI * i * ((it * factorDepth) / ((factorDepth / (it * factorDepth)))) / factorDepth));
							}
						}
						break;

					case 4:
						//std::cout << "\nno good value";
						if (graph[0] != pregrapha[0]) {
							for (int i = 0; i < graph.size(); i++) {
								pregraphs[i] += cos((M_PI * i * (factorDepth / (it * factorDepth))) / factorDepth);
								pregrapha[i] -= cos((M_PI * i * (factorDepth / (it * factorDepth))) / factorDepth);
								invpregrapha[i] -= cos((M_PI * i * ((it * factorDepth) / ((factorDepth))) / factorDepth));
								invpregraphs[i] += cos((M_PI * i * ((it * factorDepth) / ((factorDepth))) / factorDepth));
							}
						}
						break;
					}
				}
			}
		}
	}


	std::cout << "\n \n \n Done. \n \n Roots sequence nested cosine frequencies: \n";
	for (int i = 0; i < frequencies[0].size(); i++) {
		std::cout << frequencies[0][i] << ", ";
	}
	std::cout << "\n Quotient sequence nested cosine frequencies: \n";
	for (int i = 0; i < frequencies[1].size(); i++) {
		std::cout << frequencies[1][i] << ", ";
	}
	std::cout << "\n Difference sequence nested cosine frequencies: \n";
	for (int i = 0; i < frequencies[2].size(); i++) {
		std::cout << frequencies[2][i] << ", ";
	}
	std::cout << "\n Addend sequence nested cosine frequencies: \n";
	for (int i = 0; i < frequencies[3].size(); i++) {
		std::cout << frequencies[3][i] << ", ";
	}
	std::cout << "\n Product sequence nested cosine frequencies: \n";
	for (int i = 0; i < frequencies[4].size(); i++) {
		std::cout << frequencies[4][i] << ", ";
	}
	std::cout << "\n Powers sequence nested cosine frequencies: \n";
	for (int i = 0; i < frequencies[5].size(); i++) {
		std::cout << frequencies[5][i] << ", ";
	}
	std::cout << "\n \n \n Done sequencing frequency modulation approximation. Sequence of frequencies (depth " << factorDepth << "). \n\n\n"
		"Now approximating powers.";

}



bool powcheckCloser(std::vector<double>pregraph, std::vector<double>graph, std::vector<double>sequence, double factorDepth) {
	std::vector<std::vector<bool>>nearingAnswerVector(2);
	nearingAnswerVector.resize(0);
	nearingAnswerVector.resize(2);
	nearingAnswerVector[0].resize(0);
	nearingAnswerVector[1].resize(0);
	for (int im = 0; im < sequence.size(); im++) {
		if (abs(graph.at(im*factorDepth)) - (sequence[im]) < abs(pregraph.at(im * factorDepth) - (sequence[im]))) {//oh, so maybe pregraph at sequence im -1?
			for (int i = 0; i < factorDepth + 1; i++) {
				nearingAnswerVector[1].push_back(true);
			}
		}
		if (abs(graph.at(im * factorDepth)) - (sequence[im]) >= abs(pregraph.at(im * factorDepth) - (sequence[im]))) {
			for (int i = 0; i < factorDepth + 1; i++) {
				nearingAnswerVector[0].push_back(true);
			}
		}
	}

	for (int iz = 0; iz < graph.size() && iz < pregraph.size(); iz++) {
		if (round(iz / factorDepth) < sequence.size()-1) {
			if (abs(pregraph[iz] - sequence[round(iz / factorDepth)]) < abs(graph[iz] - sequence[round(iz / factorDepth)])) {
				nearingAnswerVector[0].push_back(true);
			}
			if (abs(pregraph[iz] - sequence[round(iz / factorDepth)]) >= abs(graph[iz] - sequence[round(iz / factorDepth)])) {
				nearingAnswerVector[1].push_back(true);
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
	std::vector<double>pregrapha;
	std::vector<double>pregraphs;
	std::vector<std::vector<double>>factors;


	factors.resize(0);
	factors.resize(6);

	std::cout << "\n Searching for combinations of curve power factors to fit sequence, scanning range " << cycles << " times for patterning.";
	for (int ie = 0; ie < 6; ie++) {
		std::cout << "\n " << ie << "Set processing powers, \n";
		sequence.resize(0);
		depth = seq[ie][seq[ie].size() - 1].size();
		for (int it = 0; it < seq[ie].size(); it++) {
			sequence.push_back(seq[ie][it].size());
		}
		graph.resize(0);
		pregrapha.resize(0);
		pregraphs.resize(0);
		graph.resize(factorDepth * (depth + 1));
		pregrapha.resize(factorDepth * (depth + 1));
		pregraphs.resize(factorDepth * (depth + 1));
		for(int z=0; z<cycles;z++){
			for (int it = 0; it < factorDepth * factorDepth; it++) {
				for (int i = 0; i < graph.size() && findInList(factors[ie], (it / (factorDepth))) == false; i++) {

					pregrapha[i] += pow(i / factorDepth, (it / (factorDepth)));
					pregraphs[i] -= pow(i / factorDepth, (it / (factorDepth)));
				}
				if (powcheckCloser(pregrapha, graph, sequence, factorDepth) == false && powcheckCloser(pregraphs, graph, sequence, factorDepth) == true
					&& findInList(factors[ie], (it / (factorDepth ))) == false) {
					if (factors[ie].empty() || !(factors[ie].back() <= 0)) {
						for (int i = 0; i < graph.size(); i++) {
							pregrapha[i] -= 2 * pow(i / factorDepth, (it / (factorDepth )));
							graph[i] += pow(i / factorDepth, (it / (factorDepth )));
						}
						factors[ie].push_back(-(it / (factorDepth )));
					}
					else {
						for (int i = 0; i < graph.size(); i++) {
							pregraphs[i] += pow(i / factorDepth, (it / (factorDepth * (factorDepth / 2))));
							pregrapha[i] -= pow(i / factorDepth, (it / (factorDepth * (factorDepth / 2))));
						}

					}
				}
				if (powcheckCloser(pregrapha, graph, sequence, factorDepth) == true && powcheckCloser(pregraphs, graph, sequence, factorDepth) == false
					&& findInList(factors[ie], (it / (factorDepth ))) == false) {
					if (factors[ie].empty() || !(factors[ie].back() > 0)) {
						for (int i = 0; i < graph.size(); i++) {
							pregraphs[i] += 2 * pow(i / factorDepth, (it / (factorDepth)));
							graph[i] += pow(i / factorDepth, (it / (factorDepth )));
						}
						factors[ie].push_back((it / (factorDepth )));
					}
					else {
						for (int i = 0; i < graph.size(); i++) {
							pregraphs[i] += pow(i / factorDepth, (it / (factorDepth )));
							pregrapha[i] -= pow(i / factorDepth, (it / (factorDepth)));
						}

					}


				}
				if (powcheckCloser(pregrapha, graph, sequence, factorDepth) == false && powcheckCloser(pregraphs, graph, sequence, factorDepth) == false) {
					for (int i = 0; i < graph.size(); i++) {
						pregraphs[i] += pow(i / factorDepth, (it / (factorDepth)));
						pregrapha[i] -= pow(i / factorDepth, (it / (factorDepth)));

					}
				}
				if (powcheckCloser(pregrapha, graph, sequence, factorDepth) == true && powcheckCloser(pregraphs, graph, sequence, factorDepth) == true
					&& findInList(factors[ie], (it / (factorDepth * (factorDepth / 2)))) == false) {
					for (int i = 0; i < graph.size(); i++) {
						if (powcheckCloser(pregrapha, pregraphs, sequence, factorDepth) == true && (factors[ie].empty() || !(factors[ie].back() > 0)) && findInList(factors[ie], (it / (factorDepth))) == false) {
							for (int i = 0; i < graph.size(); i++) {
								pregraphs[i] += 2 * pow(i / factorDepth, (it / (factorDepth )));
								graph[i] += pow(i / factorDepth, (it / (factorDepth)));

							}
							factors[ie].push_back((it / (factorDepth )));
						}
						if (powcheckCloser(pregraphs, pregrapha, sequence, factorDepth) == true && (factors[ie].empty() || !(factors[ie].back() <= 0)) && findInList(factors[ie], (it / (factorDepth ))) == false) {
							for (int i = 0; i < graph.size(); i++) {
								pregrapha[i] -= 2 * pow(i / factorDepth, (it / (factorDepth )));
								graph[i] -= pow(i / factorDepth, (it / (factorDepth )));

							}
							factors[ie].push_back(-(it / (factorDepth )));
						}
					}
				}
			}
		}
	}


	std::cout << "\n \n \n Done. \n \n Roots sequence curve factors: \n";
	for (int i = 0; i < factors[0].size(); i++) {
		std::cout << factors[0][i] << ", ";
	}
	std::cout << "\n Quotient sequence curve factors: \n";
	for (int i = 0; i < factors[1].size(); i++) {
		std::cout << factors[1][i] << ", ";
	}
	std::cout << "\n Difference sequence curve factors: \n";
	for (int i = 0; i < factors[2].size(); i++) {
		std::cout << factors[2][i] << ", ";
	}
	std::cout << "\n Addend sequence curve factors: \n";
	for (int i = 0; i < factors[3].size(); i++) {
		std::cout << factors[3][i] << ", ";
	}
	std::cout << "\n Product sequence curve factors: \n";
	for (int i = 0; i < factors[4].size(); i++) {
		std::cout << factors[4][i] << ", ";
	}
	std::cout << "\n Powers sequence curve factors: \n";
	for (int i = 0; i < factors[5].size(); i++) {
		std::cout << factors[5][i] << ", ";
	}

	std::cout << "\n \n \n Done sequencing powers approximation scheduling. Sequence of factors (depth " << factorDepth << "). \n";


}





int main()
{
	double depth;
	double factorizationDepth;
	bool done = false;


	while (!done) {

		std::cout << "     Welcome to UniqueSolution, a problem solver: \n" <<
			" The Number of Unique Identities Created by the Six Major Operations on Any Two Integers 0-X \n" <<
			" Enter Depth of the Six Functions sequences to be indexed: ";
		std::cin >> depth;
		std::cout << "\n Processing depth: " << depth << ".";
		std::cout << "\n And, for after the sequences are loaded:" <<
			"\n Factorization hunt depth (number of factors to consider before rerouting (suggest minimum " << depth << ")) : ";
		std::cin >> factorizationDepth;
		std::cout<<"\n Cycle length: (even numbers greater than two, suggest minimum " << round(pow(depth,.5)) << ")) : ";
		std::cin >> cycles;
		std::cout << "\n Generating sequences depth " << depth;
		calcSeq(depth);


		std::cout << "\n \n OK......  Cosine Frequencies Solutions Solver run. \n";
		modulaSeq(uniqueValues, factorizationDepth);

		std::cout << "\n \n OK......  Cosine Frequency Modulations Solutions Solver run. \n";
		fmmodulaSeq(uniqueValues, factorizationDepth);


		std::cout << "\n OK......  Relations Curve Solutions Solver run. \n";
		curveSeq(uniqueValues, factorizationDepth);

		std::cout << "\n done?(0/1): ";
		std::cin >> done;
		std::cout << "\n";
	}
	std::cin >> done;
	return !done;
}