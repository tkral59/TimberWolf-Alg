#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <random>
#include <set>
#include "Objects.hpp"
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



void perturb(std::vector<Grid*>& population, const std::map<std::string, Net>& nets, float w1, float w2, int wireConstraint) {
    std::vector<Grid*> nextGeneration;
    std::random_device rd;
    std::mt19937 gen(rd());

    for (size_t i = 0; i < population.size(); ++i) {
        bool routable = true; // Assuming this flag is needed for your cost calculation
        // Selection
        Grid* parent1 = Grid::tournamentSelection(population, 5, nets, w1, w2, wireConstraint, routable);
        Grid* parent2 = Grid::tournamentSelection(population, 5, nets, w1, w2, wireConstraint, routable);

        // Crossover - Implement your own crossover logic
        Grid* child = crossover(parent1, parent2); // You need to define how crossover works for your Grid objects

        // Mutation
        child->mutation(); // Apply mutation to the child

        // Add the new child to the next generation
        nextGeneration.push_back(new Grid(*child)); // Assuming deep copy is appropriate here
    }

    // Replace old generation with new generation
    for (auto* grid : population) {
        delete grid; // Clean up memory
    }
    population.swap(nextGeneration);
}


void createInitialGrids(const std::map<std::string, Node>& nodes, std::vector<Grid*>& gridConfigurations, int k) {

    int gridSize = std::ceil(std::sqrt(nodes.size())) * 1.1; // Adjust the scaling factor as necessary

    std::random_device rd;
    std::mt19937 g(rd());

    for (int i = 0; i < k; ++i) {
        Grid* grid = new Grid(nodes); // Using your Grid constructor that takes nodes

        // Clear existing placement if needed (assuming such a method exists)
        grid->clear(); // Make sure you implement this in Grid to reset its state

        std::vector<Node*> terminals, nonTerminals;
        for (const auto& pair : nodes) {
            if (pair.second.isTerminal())
                terminals.push_back(const_cast<Node*>(&pair.second)); // Unsafe, but needed for example. Better to adjust your design to avoid const_cast.
            else
                nonTerminals.push_back(const_cast<Node*>(&pair.second));
        }

        // Shuffle for random placement
        std::shuffle(nonTerminals.begin(), nonTerminals.end(), g);
        std::shuffle(terminals.begin(), terminals.end(), g);

        // Place terminals and non-terminals using your defined methods
        // These methods need to be defined in your Grid class to handle the placement logic
        grid->placeTerminals(terminals);
        grid->placeNonTerminals(nonTerminals);

        // Save the grid configuration
        gridConfigurations.push_back(grid);
    }
}

// Part of your finalization or analysis function
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
        
        std::lock_guard<std::mutex> lock(offspringMutex); // Protecting shared access to the offspring vector.
        offspring.push_back(child);
    }
}


void multithreadedCrossover(std::vector<Grid*>& offspring, const std::vector<Grid*>& parents, const std::map<std::string, Net>& nets, unsigned int numThreads) {
    std::vector<std::thread> threads;
    std::mutex offspringMutex; // Protects access to the offspring vector.
    
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
	 
    // Now call the multithreaded crossover instead of the single-threaded version
    multithreadedCrossover(offspring, selectedParents, nets, std::thread::hardware_concurrency());
}
*/
