#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <cmath>
#include "Objects.hpp"
#include <sstream>
#include<random>
using namespace std;

//NODE CLASS

Node::Node() {
    name = "";
    vector<Net*> netlist;
    nets = netlist;
    xcoord = 0;
    ycoord = 0;
}
Node::Node(string name, vector<Net*> nets, int xcoord, int ycoord) {
    this->name = name;
    this->nets = nets;
    this->xcoord = xcoord;
    this->ycoord = ycoord;
}

Node::~Node() {
    // Cleanup if needed
}
std::vector<Net*> Node::getNets() const {
    return nets;
}
void Node::addNet(Net* net) {
    nets.push_back(net);
}

void Node::removeNet(Net* net) {
    nets.erase(std::remove(nets.begin(), nets.end(), net), nets.end());
    for (size_t i = 0; i < nets.size(); ++i) {
        if (nets[i] == net) {
            nets.erase(nets.begin() + i);
            break;
        }
    }
}

string Node::getName() const {
    return this->name;
}

void Node::setCoords(int x, int y) {
    xcoord = x;
    ycoord = y;
}

//SQUARE CLASS

square::square() {
    type = squareType::Routing;
    node = nullptr;
    wires = 0;
}

square::square(squareType type, Node* n, int wires) {
    this->type = type;
    node = n;
    this->wires = wires;
}

void square::setType(squareType s) {
    type = s;
}

squareType square::getType() {
    return type;
}

void square::incWires() {
    wires++;
}

void square::decWires() {
    wires--;
}

void square::setNode(Node* n) {
    node = n;
}

bool square::PinsFull() {
    if (pins[0] == nullptr || pins[1] == nullptr) return false;
    else return true;
}


//UTILGRID CLASS

utilGrid::utilGrid() {

}

utilGrid::utilGrid(vector<vector<square>> ogrid) {
    int n = ogrid.size();
    std::vector<vector<square>> spacedMatrix(2 * n - 1, vector<square>(2 * n - 1));

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            spacedMatrix[2 * i][2 * j] = ogrid[i][j];
        }
    }
    for (int i = 0; i < (n*2)-1; ++i) {
        for (int j = 0; j < (n * 2) - 1; ++j) {
            if (!(spacedMatrix[i][j].getType() == squareType::Node || spacedMatrix[i][j].getType() == squareType::Terminal)) {
                square s = square(squareType::Routing, 0);
                spacedMatrix[i][j] = s;
            }
        }
    }
    grid = spacedMatrix;
}

void utilGrid::write(int x, int y, square s) {
    if (x > grid.size() || y > grid[0].size()) cout << "Write dims outside of grid." << endl;
    else grid[x][y] = s;
}

void utilGrid::swap(int x1o, int y1o, int x2o, int y2o) {
    int x1 = x1o * 2;
    int x2 = x2o * 2;
    int y1 = y1o * 2;
    int y2 = y2o * 2;
    if ((grid[x1][y1].getType() == squareType::Node && grid[x2][y2].getType() == squareType::Node) || (grid[x1][y1].getType() == squareType::Terminal && grid[x2][y2].getType() == squareType::Terminal)) { // if appropiate type and matching
        if (x1 > grid.size() || y1 > grid[0].size() || x2 > grid.size() || y2 > grid[0].size()) cout << "Swap dims outside of grid." << endl;
        else {
            square temp = grid[x1][y1];
            grid[x1][y1] = grid[x2][y2];
            grid[x2][y2] = temp;
        }
    }
    else if (grid[x1][y1].getType() == squareType::Node || grid[x2][y2].getType() == squareType::Node || grid[x1][y1].getType() == squareType::Terminal || grid[x2][y2].getType() == squareType::Terminal) cout << "Error: Trying to swap non-matching types." << endl;
    else cout << "Error: Trying to swap either Blank or Routing type square." << endl;
}

int utilGrid::calcCost() {
    // Implement the function to calculate cost
}


//GRID CLASS

// Adjust Grid constructor to automatically calculate dimensions and perform initial placement
Grid::Grid(const std::map<std::string, Node>& nodes) {
    int totalNodes = nodes.size();
    // Estimate grid size: square root of total nodes times a factor (>1) for spacing, rounded up
    int gridSize = ceil(sqrt(totalNodes) * 1.1); // Adjust the 1.1 factor as needed

    // Initialize an empty grid
    grid = vector<vector<square>>(gridSize, vector<square>(gridSize, square()));

    // Perform initial placement
    initialPlacement(nodes);
}

void Grid::write(int x, int y, square s) {
    grid[x][y] = s;
    ug.write(x * 2, y * 2, s);
}

void Grid::swap(int x1, int y1, int x2, int y2) {
    if ((grid[x1][y1].getType() == squareType::Node && grid[x2][y2].getType() == squareType::Node) || (grid[x1][y1].getType() == squareType::Terminal && grid[x2][y2].getType() == squareType::Terminal)) { // if appropiate type and matching
        if (x1 > grid.size() || y1 > grid[0].size() || x2 > grid.size() || y2 > grid[0].size()) cout << "Swap dims outside of grid." << endl;
        else {
            square temp = grid[x1][y1];
            grid[x1][y1] = grid[x2][y2];
            grid[x2][y2] = temp;
            ug.swap(x1, y1, x2, y2);
        }
    }
    else if (grid[x1][y1].getType() == squareType::Node || grid[x2][y2].getType() == squareType::Node || grid[x1][y1].getType() == squareType::Terminal || grid[x2][y2].getType() == squareType::Terminal) cout << "Error: Trying to swap non-matching types." << endl;
}

void Grid::initialPlacement(const std::map<std::string, Node>& nodes) {
    std::vector<std::pair<std::string, Node>> terminals, nonTerminals;
    // Separate terminals and non-terminals
    for (const auto& node : nodes) {
        if (node.second.isTerminal()) terminals.push_back(node);
        else nonTerminals.push_back(node);
    }

    // Shuffle terminals for random perimeter placement
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(terminals.begin(), terminals.end(), g);

    // Place terminals on the perimeter
    int perimeterCount = 0, gridSize = grid.size();
    int maxPerimeterPositions = (gridSize - 1) * 4; // Calculate available perimeter positions
    for (const auto& terminal : terminals) {
        if (perimeterCount < maxPerimeterPositions) {
            // Calculate position (x,y) based on perimeterCount
            int x = perimeterCount % (gridSize - 1);
            int y = perimeterCount / (gridSize - 1);
            if (y == 1) x = gridSize - 1;
            else if (y == 2) x = gridSize - 1 - x;
            else if (y == 3) y = gridSize - 1;
            else y = 0;

            write(x, y, square(squareType::Terminal, &terminal.second)); // Adjust for your implementation
            perimeterCount++;
        }
    }

    // Shuffle and place non-terminals within the grid's interior
    std::shuffle(nonTerminals.begin(), nonTerminals.end(), g);
    for (const auto& nonTerminal : nonTerminals) {
        bool placed = false;
        for (int i = 1; i < gridSize - 1 && !placed; i++) {
            for (int j = 1; j < gridSize - 1 && !placed; j++) {
                if (grid[i][j].getType() == squareType::Routing) { // Assuming Routing indicates an empty spot
                    write(i, j, square(squareType::Node, &nonTerminal.second)); // Adjust for your implementation
                    placed = true;
                }
            }
        }
    }
}


int Grid::calcCost() const {
    int totalCost = 0;
    for (const auto& netPair : nets) { // Assuming 'nets' is accessible and stores the Net objects
        const Net& net = netPair.second;

        // Calculate the wirelength for this net by summing the Manhattan distances
        // between all pairs of nodes connected by this net.
        for (size_t i = 0; i < net.Nodes.size(); ++i) {
            for (size_t j = i + 1; j < net.Nodes.size(); ++j) {
                int manhattanDistance = abs(net.Nodes[i]->getX() - net.Nodes[j]->getX()) +
                                        abs(net.Nodes[i]->getY() - net.Nodes[j]->getY());
                totalCost += manhattanDistance;
            }
        }
    }
    return totalCost;
}
/* Can use anyone which is better
int Grid::getCost() {
    int totalCost = 0;
    
    // Assuming 'nets' is a member of Grid that stores all the nets and their connected nodes
    for (const auto& netPair : nets) {
        const Net& net = netPair.second;
        
        // Calculate the total wirelength for this net
        for (size_t i = 0; i < net.Nodes.size(); ++i) {
            for (size_t j = i + 1; j < net.Nodes.size(); ++j) {
                int dx = abs(net.Nodes[i]->getX() - net.Nodes[j]->getX());
                int dy = abs(net.Nodes[i]->getY() - net.Nodes[j]->getY());
                totalCost += dx + dy; // Add the Manhattan distance to the total cost
            }
        }
    }
    
    return totalCost;
}
*/
