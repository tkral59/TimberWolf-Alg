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
#include <sstream> // Add this at the top of Scanner.cpp

using namespace std;

void read(const string netfile, const string nodefile, const string plfile, map<string, Node>& nodes, map<string, Net>& nets, int& numNets, int& numPins, int& numNodes, int& numTerminals) {
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
					nets[newnet.name].isCritical = false;

					netcount++;
				}
				else {
					nets["n" + to_string(netcount - 1)].Nodes.push_back(&nodes[words[0]]);
					nets["n" + to_string(netcount - 1)].nodesSize++;
					nodes[words[0]].addNet(&nets["n" + to_string(netcount - 1)]);
				}
			}
		}
		myFile.close();
	}
	else cout << "Failed to open net file." << endl;

	//read plfile
	myFile.open(plfile);
	iter = 0;
	if (myFile.is_open()) {
		while (getline(myFile, line)) {
			if (iter < 6) {
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

				nodes[words[0]].setXY(stoi(words[1]), stoi(words[2]));

			}
		}
		myFile.close();
	}
}

bool sortByValueDescending(const std::pair<string, Net>& a, const std::pair<string, Net>& b) {//ChatGPT
	return a.second.nodesSize > b.second.nodesSize;
}

void findTop10Percent(const std::map<string, Net>& inputMap) {//ChatGPT
	// Calculate the number of elements that constitute the top 10%
	int top10PercentSize = inputMap.size() * 0.1;

	// Convert the map to a vector of pairs for sorting
	std::vector<std::pair<string, Net>> vec(inputMap.begin(), inputMap.end());

	// Sort the vector by value in descending order
	std::sort(vec.begin(), vec.end(), sortByValueDescending);

	// Output the top 10% elements
	std::cout << "Top 10% elements:" << std::endl;
	for (int i = 0; i < top10PercentSize; ++i) {
		vec[i].second.isCritical = true;
	}
}
bool isNumeric(const std::string& str) {
	return !str.empty() && std::all_of(str.begin(), str.end(), [](char c) { return std::isdigit(c) || c == '-'; });
}


vector<Result> createInitialGrids(const std::map<std::string, Node>& nodes, int k, float const w1, float const w2, float const w3, map<string, Net> const nets, int wireConstraint) {
	vector<Result> init;
	std::random_device rd;
	std::mt19937 g(rd());

	for (int i = 0; i < k; ++i) {
		Grid grid(nodes); // Use `grid` instead of `g` to avoid confusion with `std::mt19937 g`
		bool routable = false;
		vector<Bounds> bounds;
		float cost = grid.calcCost(w1, w2, w3, nets, routable, wireConstraint, bounds);

		init.emplace_back(std::move(grid), cost, routable, std::move(bounds));
	}
	cout << "Finished creating Initial Grids" << endl;
	return init;
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

void exportForVisualization(Result r, const std::string& filename) {
	std::ofstream file(filename);
	file << "Net,Xmin,Ymin,Xmax,Ymax\n";
	for (auto b : r.bounds) {
		// Assume bounds are calculated and stored somewhere accessible
		file << b.net->name << "," << b.x1 << "," << b.y1 << "," << b.x2 << "," << b.y2 << "\n";
	}
	file.close();
}
/*
void performCrossoversThread(std::vector<Grid*>& offspring, const std::vector<Grid*>& parents, const std::map<std::string, Net>& nets, int startIdx, int endIdx, std::mutex& offspringMutex, map<string, Node> nodes) {
	for (int i = startIdx; i < endIdx && (i + 1) < parents.size(); i += 2) {
		// Perform crossover on parents[i] and parents[i+1]
		Grid* child = crossover(parents[i], parents[i + 1], nets, nodes); // Ensure your crossover function is thread-safe.
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
*/
bool compareByFloat(const Result& a, const Result& b) { //ChatGPT
	return a.cost > b.cost; // Change to < for ascending order
}

vector<Result> tournamentSelection(std::vector<Result> population, const std::map<std::string, Net>& nets, float percentage) { //ChatGPT
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> dist(0, population.size() - 1);

	int topKNum = population.size() * percentage;
	vector<Result> topFive = population;
	sort(topFive.begin(), topFive.end(), compareByFloat);
	auto h = topFive.begin();
	advance(h, topKNum);
	vector<Result> newTopFive(topFive.begin(), h);
	return newTopFive;
}

vector<Result> perturb(std::vector<Result>& population, const std::map<std::string, Net>& nets, float w1, float w2, float w3, int wireConstraint, map<string, Node> nodes, float selectP, float crossP, float mutP) {
	std::vector<Result> nextGeneration;
	std::random_device rd;
	std::mt19937 gen(rd());

	int count = 0; //to keep track of how many offspring created already

	//Selection:
	nextGeneration = tournamentSelection(population, nets, selectP);
	count += nextGeneration.size();
	//Crossover
	int cCount = population.size() * crossP;
	for (int i = 0; i < cCount; ++i) {
		int i1 = rand() % count;
		int i2 = rand() % count;
		Grid* parent1 = &nextGeneration[i1].g;
		Grid* parent2 = &nextGeneration[i2].g;
		Grid* child = crossover(parent1, parent2, nets, nodes); //may be creation error here
		bool rout = false;
		vector<Bounds> b;
		float cost = child->calcCost(w1, w2, w3, nets, rout, wireConstraint, b);
		Result r(*child, cost, rout, b);
		nextGeneration.push_back(r);

	}
	count += cCount;
	//Mutation
	for (int i = count - 1; i < population.size(); i++) {
		int in = rand() % count;
		Grid* copy = &nextGeneration[in].g;
		int rx = rand() % copy->getGridX();
		int ry = rand() % copy->getGridY();
		copy->mutation(rx, ry);
		bool rout = false;
		vector<Bounds> b;
		float cost = copy->calcCost(w1, w2, w3, nets, rout, wireConstraint, b);
		Result r(*copy, cost, rout, b);
		nextGeneration.push_back(r);
	}
	return nextGeneration;
}

double generateInitialTemp(vector<Result> init, double prob, float const w1, float const w2, float const w3, map<string, Net> const nets, bool& routable, int wireConstraint) {
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
				Grid copy = r.g;
				int rx = rand() % copy.getGridX(); //create initial grid x param;
				int ry = rand() % copy.getGridY();//create initial grid y param;
				copy.mutation(rx, ry);
				bool route = false;
				vector<Bounds> b;
				float cost = copy.calcCost(w1, w2, w3, nets, routable, wireConstraint, b);
				Result n(copy, cost, route, b);
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

Result simulatedAnnealing(vector<Result> initialGrids, float const w1, float const w2, float const w3, map<string, Net> const nets, int wireConstraint, map<string, Node> nodes) {
	bool routable = false; // Ensure it's declared
	double t = generateInitialTemp(initialGrids, 5., w1, w2, w3, nets, routable, wireConstraint);
	double initT = t;
	vector<Result> population = initialGrids;
	vector<Result> new_pop;
	vector<Result> best_pop;
	double deltaC = 0;
	cout << "Initial Cost: " << bestCost(population).cost << endl;
	int iteration = 1;
	while (t > 0) {
		while (routable == false) {
			cout << "Iteration " << iteration << "; Temp = " << t << endl;
			new_pop = perturb(population, nets, w1, w2, w3, wireConstraint, nodes, 0.5, 0.25, 0.25); //NEED PERTURB FUNCTION //NEEDS TO RETURN LIST OF GRIDS : COST : ROUTABLE?
			Result nbc = bestCost(new_pop);
			deltaC = nbc.cost - bestCost(population).cost; //NEED BEST COST FUNCTION
			cout << "\t Delta C = " << deltaC << endl;

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
				if (nbc.cost < bestCost(best_pop).cost) {
					best_pop = new_pop;
				}
				cout << "\t \t new best population!" << endl;
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
void simulatedAnnealing(Grid& initialGrid, const std::map<std::string, Net>& nets, int wireConstraint, float initialTemperature = 100.0f, int totalSteps = 10000) {
	srand(static_cast<unsigned>(time(nullptr))); // Seed the RNG

	float temperature = initialTemperature;
	Grid currentGrid = initialGrid;
	bool routable;
	std::vector<Bounds> bounds;
	float currentCost = currentGrid.calcCost(1.0, 1.0, 1.0, nets, routable, wireConstraint, bounds);

	Grid bestGrid = currentGrid;
	float bestCost = currentCost;

	for (int step = 0; step < totalSteps && temperature > 1e-3; ++step) {
		Grid newGrid = currentGrid; // Create a new candidate by copying the current grid
		newGrid.mutation(rand() % newGrid.getGridX(), rand() % newGrid.getGridY()); // Apply mutation

		float newCost = newGrid.calcCost(1.0, 1.0, 1.0, nets, routable, wireConstraint, bounds); // Calculate the new cost

		if ((newCost < currentCost) || (exp((currentCost - newCost) / temperature) > static_cast<float>(rand()) / RAND_MAX)) {
			currentGrid = newGrid;
			currentCost = newCost;

			if (newCost < bestCost) {
				bestGrid = newGrid;
				bestCost = newCost;
			}
		}

		// Update temperature based on a simple cooling schedule
		temperature *= 0.95;
	}

	// bestGrid now contains the optimized grid configuration
}
*/

int main() {
	std::string netfile = "P2Benchmarks\\ibm01\\ibm01.nets";
	std::string nodefile = "P2Benchmarks\\ibm01\\ibm01.nodes";
	std::string plfile = "P2Benchmarks\\ibm01\\ibm01.pl";

	std::map<std::string, Node> nodes;
	std::map<std::string, Net> nets;
	int numNet = 0, numPins = 0, numNode = 0, numTerminals = 0;

	try {
		read(netfile, nodefile, plfile, nodes, nets, numNet, numPins, numNode, numTerminals);
		std::vector<Result> init = createInitialGrids(nodes, 10, 1.0, 1.0, 1.0, nets, 4);
		Result bestResult = simulatedAnnealing(init, 1.0, 1.0, 1.0, nets, 4, nodes);

		std::cout << "Optimization completed. Displaying results:\n";
		// Assuming you now return a single best result from simulatedAnnealing
		std::cout << "Best cost: " << bestResult.cost << std::endl;
	}
	catch (const std::exception& e) {
		std::cerr << "Error during optimization: " << e.what() << std::endl;
	}

	return 0;
}