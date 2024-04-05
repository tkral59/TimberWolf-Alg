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

void readPLFileAndUpdateNodes(const std::string& filename, std::map<std::string, Node>& nodes) {
    std::ifstream plFile(filename);

    // Directly check if the file is open and readable
    if (!plFile.is_open()) {
        std::cerr << "Error: .pl file '" << filename << "' does not exist or cannot be opened." << std::endl;
        return;  // Stop processing if the file couldn't be opened
    }

    std::string line;
    while (getline(plFile, line)) {
        if (line.empty() || line[0] == '#' || line.find("UCLA") != std::string::npos) continue; // Skip comments and headers

        std::istringstream iss(line);
        std::string nodeName, temp;
        int x, y;
        if (iss >> nodeName >> x >> y >> temp) { // Assuming 'temp' captures any trailing attributes, which are ignored here
            // Check if nodeName exists in the nodes map and update positions
            auto it = nodes.find(nodeName);
            if (it != nodes.end()) {
                it->second.x = x;
                it->second.y = y;
                // Optionally, mark node as terminal based on additional logic
                it->second.isTerminal = true; // Example, adjust as needed
            } else {
                // If the nodeName does not exist in the map, you might want to add it,
                // Or handle it depending on your application's needs
                // nodes[nodeName] = Node(nodeName, x, y, true); // Example, adjust as needed
            }
        }
    }

    // Confirm successful reading:
    std::cout << "Successfully read and processed '" << filename << "'." << std::endl;
}


vector<Result> createInitialGrids(const std::map<std::string, Node>& nodes, int k, float const w1, float const w2, map<string, Net> const nets, int wireConstraint) {

    vector<Result> init;
    std::random_device rd;
    std::mt19937 g(rd());

    for (int i = 0; i < k; ++i) {
        Grid grid(nodes); // Directly construct Grid object
        grid.initialPlacement(nodes);
        bool routable = false;
        vector<Bounds> bounds;
        float cost = grid.calcCost(w1, w2, nets, routable, wireConstraint, bounds);
        
        // Use std::move to move the Grid object into the Result, avoiding unnecessary copies
        init.emplace_back(std::move(grid), cost, routable, std::move(bounds));
    }

    return init; // Don't forget to return the populated vector
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

void performCrossoversThread(std::vector<Grid*>& offspring, const std::vector<Grid*>& parents, const std::map<std::string, Net>& nets, int startIdx, int endIdx, std::mutex& offspringMutex, map<string, Node> nodes) {
	for (int i = startIdx; i < endIdx && (i + 1) < parents.size(); i += 2) {
		// Perform crossover on parents[i] and parents[i+1]
		Grid* child = crossover(parents[i], parents[i + 1], nets, nodes); // Ensure your crossover function is thread-safe.

		std::lock_guard<mutex> lock(offspringMutex); // Protecting shared access to the offspring vector.
		offspring.push_back(child);
	}
}


void multithreadedCrossover(std::vector<Grid*>& offspring, const std::vector<Grid*>& parents, const std::map<std::string, Net>& nets, unsigned int numThreads, std::map<string, Node> nodes) {
    vector<thread> threads;
    mutex offspringMutex; // Protects access to the offspring vector.

    int segmentSize = parents.size() / numThreads; // Determine workload size per thread.

    for (unsigned int i = 0; i < numThreads; ++i) {
        int startIdx = i * segmentSize;
        int endIdx = (i == numThreads - 1) ? parents.size() : (i + 1) * segmentSize; // Ensure the last thread covers any remaining parents.

        // Launch a thread to process its segment of the parents vector using a lambda to capture the necessary variables.
        threads.emplace_back([&, startIdx, endIdx, i]() {
            performCrossoversThread(offspring, parents, nets, startIdx, endIdx, offspringMutex, nodes);
        });
    }

    // Wait for all threads to complete their tasks.
    for (std::thread& t : threads) {
        if (t.joinable()) {
            t.join();
        }
    }
}


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

vector<Result> perturb(std::vector<Result>& population, const std::map<std::string, Net>& nets, float w1, float w2, int wireConstraint, map<string, Node> nodes, float selectP, float crossP, float mutP) {
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
		float cost = child->calcCost(w1, w2, nets, rout, wireConstraint, b);
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
		float cost = copy->calcCost(w1, w2, nets, rout, wireConstraint, b);
		Result r(*copy, cost, rout, b);
		nextGeneration.push_back(r);
	}
	return nextGeneration;
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
				Grid copy = r.g;
				int rx = rand() % copy.getGridX(); //create initial grid x param;
				int ry = rand() % copy.getGridY();//create initial grid y param;
				copy.mutation(rx, ry);
				bool route;
				vector<Bounds> b;
				float cost = copy.calcCost(w1, w2, nets, routable, wireConstraint, b);
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

Result simulatedAnnealing(vector<Result> initialGrids, float const w1, float const w2, map<string, Net> const nets, int wireConstraint, map<string, Node> nodes) {
	bool routable = false;
	double t = generateInitialTemp(initialGrids, 5., w1, w2, nets, routable, wireConstraint); //initial temp
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
			new_pop = perturb(population, nets, w1, w2, wireConstraint, nodes, 0.5, 0.25, 0.25); //NEED PERTURB FUNCTION //NEEDS TO RETURN LIST OF GRIDS : COST : ROUTABLE?
			deltaC = bestCost(new_pop).cost - bestCost(population).cost; //NEED BEST COST FUNCTION
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
				best_pop = new_pop;
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

void simulatedAnnealing(Grid& initialGrid, const std::map<std::string, Net>& nets, int wireConstraint, float initialTemperature = 100.0f, int totalSteps = 10000) {
    srand(static_cast<unsigned>(time(nullptr))); // Seed the RNG

    float temperature = initialTemperature;
    Grid currentGrid = initialGrid;
    bool routable;
    std::vector<Bounds> bounds;
    float currentCost = currentGrid.calcCost(1.0, 1.0, nets, routable, wireConstraint, bounds);

    Grid bestGrid = currentGrid;
    float bestCost = currentCost;

    for (int step = 0; step < totalSteps && temperature > 1e-3; ++step) {
        Grid newGrid = currentGrid; // Create a new candidate by copying the current grid
        newGrid.mutation(rand() % newGrid.getGridX(), rand() % newGrid.getGridY()); // Apply mutation

        float newCost = newGrid.calcCost(1.0, 1.0, nets, routable, wireConstraint, bounds); // Calculate the new cost

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

int main() {
    // Paths to your net and node files need to be specified here
    std::string netfile = "/Users/Karan/Downloads/vlsi3/ibmISPD02Bench_Bookshelf/ibm01/ibm01.nets"; // Update with the actual path to your net file
    std::string nodefile = "/Users/Karan/Downloads/vlsi3/ibmISPD02Bench_Bookshelf/ibm01/ibm01.nodes";
	std::string filename = "path/to/your/file.pl";
    std::map<std::string, Node> nodes;

    // Attempt to read the .pl file and update node positions
    readPLFileAndUpdateNodes(filename, nodes);
 // Update with the actual path to your node file

    // Containers for nodes and nets populated by the read function
    std::map<std::string, Node> nodes;
    std::map<std::string, Net> nets;
    int numNet = 0, numPins = 0, numNode = 0, numTerminals = 0;

    // Reading node and net data from files
    read(netfile, nodefile, nodes, nets, numNet, numPins, numNode, numTerminals);

    // Generating initial grid configurations for optimization
    std::vector<Result> init = createInitialGrids(nodes, 10, 0.5, 0.5, nets, 4);

    // Executing the simulated annealing optimization process
    simulatedAnnealing(init, 0.5, 0.5, nets, 4, nodes);

    // Following the optimization, you might want to output the optimized results
    // This section is left generic since the specifics depend on your Result structure and desired outputs
    std::cout << "Optimization completed. Displaying results:" << std::endl;
    // Example: Displaying the cost of the best result
    if (!init.empty()) {
        std::cout << "Best cost: " << init.front().cost << std::endl;
    } else {
        std::cout << "No results available." << std::endl;
    }

    return 0; // Ensure to return an integer type for main()
}
