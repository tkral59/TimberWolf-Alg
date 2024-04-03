#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <random>
#include <set>
#include <cmath>
#include <thread>
#include <mutex>
#include "Objects.hpp"
//#include "Scanner.hpp"
#include <sstream> // Add this at the top of Scanner.cpp

using namespace std;

void read(const string netfile, const string nodefile, map<string, Node>& nodes, map<string, Net>& nets, int& numNets, int& numPins, int& numNodes, int& numTerminals) {
	nodes.clear();
	nets.clear();

	string line;
	fstream myFile;

	//NODE FILE
	int iter = 0;
	myFile.open(nodefile);
	if (myFile.is_open()) {
		while (getline(myFile, line)) {
			if (iter < 5) iter++;
			else {
				vector<string> words;
				string temp;
				stringstream ss(line);

				while (ss >> temp) {
					words.push_back(temp); //delimits string into words by spaces
				}

				if (words[0] == "NumNodes") numNodes = stoi(words[2]);
				else if (words[0] == "NumTerminals") numTerminals = stoi(words[2]);
				else {
					vector<Net*> newlist;
					Node newnode(words[0], newlist);
					nodes[words[0]] = newnode;
				}
			}
		}
		myFile.close();
	}
	else {
		cout << "Failed to open node file." << endl;
	}


	//NETFILE
	myFile.open(netfile);
	iter = 0;
	int netcount = 0;
	if (myFile.is_open()) {
		while (getline(myFile, line)) {
			if (iter < 5) {
				iter++;
				continue;
			}
			else {
				vector<string> words;
				string temp;
				stringstream ss(line);

				while (ss >> temp) {
					words.push_back(temp); //delimits string into words by spaces
				}

				if (words[0] == "NumNets") numNets = stoi(words[2]);
				else if (words[0] == "NumPins") numPins = stoi(words[2]);
				else if (words[0] == "NetDegree") {
					Net newnet("n" + to_string(netcount));
					nets[newnet.name] = newnet;
					netcount++;
				}
				else {
					nets["n" + to_string(netcount - 1)].Nodes.push_back(&nodes[words[0]]);
					nodes[words[0]].addNet(&nets["n" + to_string(netcount - 1)]);
				}
			}
		}
	}
	else cout << "Failed to open net file." << endl;

}

vector<Result> createInitialGrids(const std::map<std::string, Node>& nodes, int k, float const w1, float const w2, map<string, Net> const nets, int wireConstraint) {

	int gridSize = std::ceil(std::sqrt(nodes.size())) * 1.1; // Adjust the scaling factor as necessary
	vector<Result> init;
	std::random_device rd;
	std::mt19937 g(rd());

	for (int i = 0; i < k; ++i) {
		//Grid* grid = new Grid(nodes); // Using your Grid constructor that takes nodes
		Grid g(nodes);
		g.initialPlacement(nodes);
		bool routable = false;
		float cost = g.calcCost(w1, w2, nets, routable, wireConstraint);
		Result r(g, cost, routable);

		// Save the grid configuration
		init.push_back(r);
	}
}

Result bestCost(vector<Result> results) {
	Result best;
	double mincost = results.at(0).cost;
	if (results.at(0).routable) return results.at(0);
	for (int i = 1; i < results.size(); i++) {
		if (results.at(i).routable) {
			return results.at(i);
		}
		else if (results.at(i).cost < mincost) {
			mincost = results.at(i).cost;
			best = results.at(i);
		}
	}
	return best;
}

Grid* crossover(Grid* parent1, Grid* parent2, const std::map<std::string, Net>& nets, map<string, Node> nodes) {
	// Assume Grid has a constructor that takes the size and nets to initialize an empty grid.
	auto child = new Grid(nodes);

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> dis(0, parent1->getGridSize() - 1);

	int crossoverPoint = dis(gen);

	// Copy up to the crossover point from parent1
	for (int i = 0; i <= crossoverPoint; ++i) {
		for (int j = 0; j < parent1->getGridSize(); ++j) {
			auto node = parent1->getSquare(i, j).getNode();
			if (node && !child->isNodePlaced(node)) {
				child->placeNode(i, j, node);
			}
		}
	}

	// Fill in the rest from parent2, avoiding duplicates
	for (int i = crossoverPoint + 1; i < parent2->getGridSize(); ++i) {
		for (int j = 0; j < parent2->getGridSize(); ++j) {
			auto node = parent2->getSquare(i, j).getNode();
			if (node && !child->isNodePlaced(node)) {
				child->placeNode(i, j, node);
			}
		}
	}

	return child;
}

void exportForVisualization(const std::map<std::string, Net>& nets, const std::string& filename) {
	std::ofstream file(filename);
	file << "Net,Xmin,Ymin,Xmax,Ymax\n";
	for (const auto& netPair : nets) {
		const Net& net = netPair.second;
		// Assume bounds are calculated and stored somewhere accessible
		file << net.name << "," << net.xmin << "," << net.ymin << "," << net.xmax << "," << net.ymax << "\n";
	}
	file.close();
}

void performCrossoversThread(std::vector<Grid*>& offspring, const std::vector<Grid*>& parents, const std::map<std::string, Net>& nets, int startIdx, int endIdx, std::mutex& offspringMutex) {
	for (int i = startIdx; i < endIdx && (i + 1) < parents.size(); i += 2) {
		// Perform crossover on parents[i] and parents[i+1]
		Grid* child = crossover(parents[i], parents[i + 1], nets); // Ensure your crossover function is thread-safe.

		std::lock_guard<mutex> lock(offspringMutex); // Protecting shared access to the offspring vector.
		offspring.push_back(child);
	}
}


void multithreadedCrossover(std::vector<Grid*>& offspring, const std::vector<Grid*>& parents, const std::map<std::string, Net>& nets, unsigned int numThreads) {
	vector<thread> threads;
	mutex offspringMutex; // Protects access to the offspring vector.

	int segmentSize = parents.size() / numThreads; // Determine workload size per thread.

	for (unsigned int i = 0; i < numThreads; ++i) {
		int startIdx = i * segmentSize;
		int endIdx = (i == numThreads - 1) ? parents.size() : (i + 1) * segmentSize; // Ensure the last thread covers any remaining parents.

		// Launch a thread to process its segment of the parents vector.
		threads.emplace_back(performCrossoversThread, std::ref(offspring), std::ref(parents), std::ref(nets), startIdx, endIdx, std::ref(offspringMutex));
	}

	// Wait for all threads to complete their tasks.
	for (std::thread& t : threads) {
		if (t.joinable()) {
			t.join();
		}
	}
}


void perturb(std::vector<Grid*>& population, const std::map<std::string, Net>& nets, float w1, float w2, int wireConstraint, map<string, Node> nodes) {
	std::vector<Grid*> nextGeneration;
	std::random_device rd;
	std::mt19937 gen(rd());

	for (size_t i = 0; i < population.size(); ++i) {
		bool routable = true; // Assuming this flag is needed for your cost calculation
		// Selection
		Grid* parent1 = Grid::tournamentSelection(population, 5, nets, w1, w2, wireConstraint, routable);
		Grid* parent2 = Grid::tournamentSelection(population, 5, nets, w1, w2, wireConstraint, routable);

		// Crossover - Implement your own crossover logic
		Grid* child = crossover(parent1, parent2, nets, nodes); // You need to define how crossover works for your Grid objects

		// Mutation
		child->mutation(0, 0); // Apply mutation to the child

		// Add the new child to the next generation
		nextGeneration.push_back(new Grid(*child)); // Assuming deep copy is appropriate here
	}

	// Replace old generation with new generation
	for (auto* grid : population) {
		delete grid; // Clean up memory
	}
	population.swap(nextGeneration);
}

double generateInitialTemp(vector<Result> init, double prob, float const w1, float const w2, map<string, Net> const nets, bool& routable, int wireConstraint) {
	double emax = 0., emin = 0.;
	double esum = 0;
	vector<double> es;
	for (Result r : init) {
		double e = 0.;
		while (e <= 0.) { //get an positive transition
			int ra = rand() % 3;
			int ri = rand() % init.size();
			if (ra <= 1) { //select
				e = r.cost - init.at(ri).cost;
			}
			else if (ra == 2) {
				//crossover
			}
			else {
				int rx = rand()%; //create initial grid x param;
				int ry = rand()%; //create initial grid y param;
				Grid copy = r.g;
				copy.mutation(rx, ry);
				bool route;
				Result n(copy, copy.calcCost(w1, w2, nets, routable, wireConstraint), route);
				e = n.cost - r.cost;
			}
		}
		es.push_back(e);
		esum += e;
	}
	auto emaxit = max_element(es.begin(), es.end());
	auto eminit = min_element(es.begin(), es.end());
	emax = *emaxit;
	emin = *eminit;
	double t = emax;
	double xo = 0.8;
	double x = 0.;
	while (x < xo) {
		x = exp(-(emax / t)) / exp(-(emin / t));
		t = t * pow((log(x) / log(xo)), (1. / prob));
	}
	return t;
}

double schedule(double temp, double initialTemp) {
	double percentComplete = (initialTemp - temp) / initialTemp;
	if (percentComplete < 0.8 || percentComplete > 0.92) {
		return 0.95 * temp;
	}
	else return 0.8 * temp;
}

Result simulatedAnnealing(vector<Result> initialGrids, float const w1, float const w2, map<string, Net> const nets, int wireConstraint) {
	bool routable = false;
	double t = generateInitialTemp(initialGrids, 5., w1, w2, nets, routable, wireConstraint); //initial temp
	double initT = t;
	vector<Result> population = initialGrids;
	vector<Result> new_pop;
	vector<Result> best_pop;
	double deltaC = 0;
	while (t > 0) {
		while (routable == false) {
			new_pop = perturb(population); //NEED PERTURB FUNCTION //NEEDS TO RETURN LIST OF GRIDS : COST : ROUTABLE?
			deltaC = bestCost(new_pop).cost - bestCost(population).cost; //NEED BEST COST FUNCTION

			// for exploration
			random_device rd;
			mt19937 gen(rd()); //seed;
			uniform_real_distribution<double> dis(0.0, 1.0);
			double r = dis(gen);
			double e = exp(deltaC / t);
			//end exploration parameters

			//if better cost, exploit
			if (deltaC < 0) {
				population = new_pop;
				best_pop = new_pop;
			}

			//chance to explore
			else if (r > e) {
				population = new_pop;
			}
			t = schedule(t, initT);
		}
	}
	return bestCost(best_pop);
}

/*
void main() {
	string netfile = "P2Benchmarks\\ibm01\\ibm01.nets";
	string nodefile = "P2Benchmarks\\ibm01\\ibm01.nodes";
	int NumTerminals = 0, NumNodes = 0, NumNets = 0, NumPins = 0;
	map<string, Node> nodes;
	map<string, Net> nets;
	read(netfile, nodefile, nodes, nets, NumNets, NumPins, NumNodes, NumTerminals);
	for (auto pair : nets) {
		cout << pair.first << ":" << endl;
		for (auto node : pair.second.Nodes) {
			cout << "\t" << node->getName() << endl;
		}
	}
}
*/

void main() {
	string netfile = "";
	string nodefile = "";
	map<string, Node> nodes;
	map<string, Net> nets;
	int numNet, numPins, numNode, numTerminals;

	read(netfile, nodefile, nodes, nets, numNet, numPins, numNode, numTerminals);
	//simulatedAnnealing();
}
