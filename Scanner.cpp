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

void runEvolutionaryAlgorithm(std::vector<Grid*>& initialPopulation, const std::map<std::string, Net>& nets, float w1, float w2, int wireConstraint) {
    std::vector<Grid*> population = initialPopulation;
    bool routable = true;
    const size_t populationSize = population.size();
    const size_t tournamentSize = 5; // Adjust based on your scenario

    for (int generation = 0; generation < 100; ++generation) { // Example: 100 generations
        std::vector<Grid*> newGeneration;
        while (newGeneration.size() < populationSize) {
            Grid* parent = Grid::tournamentSelection(population, tournamentSize, nets, w1, w2, wireConstraint, routable);
            // Here, clone or mutate `parent` to create a new grid instance
            Grid* offspring = new Grid(*parent); // Simplified: directly cloning
            // Consider applying mutations here
            newGeneration.push_back(offspring);
        }

        // Cleanup previous generation
        for (auto* grid : population) delete grid;
        population = std::move(newGeneration);
    }

    // Cleanup final generation
    for (auto* grid : population) delete grid;
}

void createInitialGrids(map<string, Node> nodes, map<string, Net> nets) {
	//create k random grid setups to act as first generation for timberwolf algorithm.
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
