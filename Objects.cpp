#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <random>
#include <map>
#include <cmath>
#include "Objects.hpp"
#include <sstream>
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

square::square(squareType type, const Node* n, int wires) {
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
    for (int i = 0; i < (n * 2) - 1; ++i) {
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


//GRID CLASS

Grid::Grid(int rows, int cols) {
    grid = vector<vector<square>>(rows, vector<square>(cols));
    ug = utilGrid(grid);
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
    std::vector<Node> nodeVec;
    for (const auto& pair : nodes) {
        nodeVec.push_back(pair.second);
    }

    // Shuffle the vector to randomize node order
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 g(rd());   // Standard mersenne_twister_engine seeded with rd()
    std::shuffle(nodeVec.begin(), nodeVec.end(), g);

    int x = 0, y = 0; // Starting coordinates for placement.
    for (const auto& node : nodeVec) {
        // Create a square with this node. Adjust according to your square constructor.
        square s(squareType::Node, &node); // Assuming your constructor can take a pointer to Node
        this->write(x, y, s); // Place the node on the grid.

        // Update coordinates for next placement, with a simple pattern for spacing.
        x += 2; // Assuming a grid spacing pattern. Adjust as needed.
        if (x >= this->grid.size()) { // Move to the next row if we've reached the end of the current row.
            x = 0;
            y += 2;
        }

        // Check for grid bounds to ensure we don't exceed grid dimensions.
        if (y >= this->grid[0].size()) {
            std::cout << "Warning: Grid is full, not all nodes have been placed." << std::endl;
            break; // Exit if we run out of grid space.
        }
    }
}

int Grid::calcCost(float const w1, float const w2, map<string, Net> const nets, bool& routable, int wireConstraint) const {
    float totalCost = 0, totalLength = 0, overlapCount = 0;
    vector<Bounds> bounded;
    for (const auto& netPair : nets) { // Assuming 'nets' is accessible and stores the Net objects
        const Net* net = &netPair.second;
        int xmin = net->Nodes.at(0)->getX(), xmax = xmin, ymin = net->Nodes.at(0)->getX(), ymax = ymin;
        Bounds newBounds;
        // Calculate the wirelength for this net by finding the x and y bounds (half-param measure)

        for (size_t i = 0; i < net->Nodes.size(); ++i) { //iterate through nodes and find min/max x/y
            int x = net->Nodes.at(i)->getX();
            int y = net->Nodes.at(i)->getY();
            if (x < xmin) xmin = x;
            if (x > xmax) xmax = x;
            if (y < ymin) ymin = y;
            if (y > ymax) ymax = y;
            totalLength = abs(xmax - xmin) + abs(ymax - ymin);
            newBounds.x1 = xmin, newBounds.x2 = xmax, newBounds.y1 = ymin, newBounds.y2 = ymax, newBounds.net = net;
        }
        bounded.push_back(newBounds);
        int olcount = 0;
        for (Bounds const bounds : bounded) {
            if (bounds.net->name != netPair.first) {
                //overlap cost (total nets overlap)
                float x_overlap = max(0, min(newBounds.x2, bounds.x2) - max(newBounds.x1, bounds.x1)); //x overlap
                float y_overlap = max(0, min(newBounds.y2, bounds.y2) - max(newBounds.y1, bounds.y1));//y overlap
                float interArea = x_overlap * y_overlap; //total intersection area
                float area1 = abs(newBounds.x2 - newBounds.x1) * abs(newBounds.y2 - newBounds.y1); //calculating area of current net box
                float area2 = abs(bounds.x2 - bounds.x1) * abs(bounds.y2 - bounds.y1); //net j box area
                if ((interArea / min(area1, area2)) <= 0.25) { //if overlap is greater than 25%
                    olcount++;
                    overlapCount++;
                    if (olcount > wireConstraint) routable = false; //if overlapping net boxes > constraint then the design is not routable
                }
            }
        }
    }
    delete& bounded;
    float ocnorm = overlapCount / nets.size(); //normalized overlap count cost => total count of nets is max, min is zero
    float den = (ug.grid.size() * ug.grid[0].size());
    float tlnorm = totalLength / den; //normallized total length cost => total grid area * net count is max, min is ~ 1
    totalCost = (w1 * tlnorm) + (w2 * ocnorm);
    return totalCost;
}
