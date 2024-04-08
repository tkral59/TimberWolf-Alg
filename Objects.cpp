#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <cmath>
#include "Objects.hpp"
#include <sstream>
#include <random>
using namespace std;

//NODE CLASS

Node::Node() {
    name = "";
    vector<Net*> netlist;
    nets = netlist;
    x = 0;
    y = 0;
    z = 0;
}

Node::Node(string name, vector<Net*> nets, int x, int y, bool isTerminal) {
    this->name = name;
    this->nets = nets;
    this->x = x;  // Corrected from xcoord
    this->y = y;  // Corrected from ycoord
    this->isTerminalFlag = isTerminal;  // Assuming you have a member variable isTerminalFlag
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

int Node::getX() const {
    return x; // Assuming 'x' is an integer member variable of Node
}

int Node::getY() const {
    return y; // Assuming 'y' is an integer member variable of Node
}

int Node::getZ() const {
    return z;
}
void Node::setXY(int newX, int newY) {
    this->x = newX;
    this->y = newY;
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
    if (this->name.empty() || x < -33 || y < -33) {
        cout << "uh oh" << endl;
        return "uh oh";
    }
    else {
        return this->name;
    }

}


bool Node::isTerminal() const {
    return name[0] == 'p';
}

void Node::setTerminal(bool terminal) {
    this->isTerminalFlag = terminal;
}


//SQUARE CLASS

square::square() {
    type = squareType::Node;
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

void square::setNode(const Node* n) {
    node = n;
}

const Node* square::getNode() {
    return node;
}

bool square::isEmpty() {
    if (node == nullptr) return true;
    else return false;
}


//UTILGRID CLASS

utilGrid::utilGrid() {

}

utilGrid::utilGrid(vector<vector<square>> ogrid) {//chatGPT aided
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

void utilGrid::move(int x1o, int y1o, int x2o, int y2o) {
    int x1 = x1o * 2;
    int x2 = x2o * 2;
    int y1 = y1o * 2;
    int y2 = y2o * 2;

    if ((grid[x1][y1].getType() == squareType::Node && grid[x2][y2].getType() == squareType::Node) || (grid[x1][y1].getType() == squareType::Terminal && grid[x2][y2].getType() == squareType::Terminal)) { // if one is node and the other is empty
        if (x1 > grid.size() || y1 > grid[0].size() || x2 > grid.size() || y2 > grid[0].size()) cout << "movement dims outside of grid." << endl;
        else {
            if (grid[x2][y2].getNode() == nullptr) { grid[x2][y2] = grid[x1][y1]; }
            else if (grid[x1][y1].getNode() == nullptr) { grid[x1][y1] = grid[x2][y2]; }
        }
    }
    else if (grid[x1][y1].getType() == squareType::Node || grid[x2][y2].getType() == squareType::Node || grid[x1][y1].getType() == squareType::Terminal || grid[x2][y2].getType() == squareType::Terminal) cout << "Error: Trying to move non-matching types." << endl;
}

//GRID CLASS

Grid::Grid() {
    // Initialization or other actions
}
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

Grid::Grid(int totalNodes) {
    int gridSize = ceil(sqrt(totalNodes) * 1.1);
    grid = vector<vector<square>>(gridSize, vector<square>(gridSize, square()));
}

void Grid::write(int x, int y, square s) {
    if (x >= 0 && x < grid.size() && y >= 0 && y < grid[x].size()) {
        grid[x][y] = s;
        Coords c(x, y);
        nodeCoords[s.getNode()->getName()] = c;
        if (2 * x < ug.grid.size() && 2 * y < ug.grid[0].size()) { // Assuming ug.grid is public; adjust if it's private
            ug.write(x * 2, y * 2, s);
        }
    }
    else {
        // Handle out-of-bounds access appropriately
    }
}

void Grid::updateEmpties(int x1, int y1, int x2, int y2, bool isTerminal) {
    Coords a(x2, y2);
    Coords b(x1, y1);
    enodes.push_back(b);
    grid[x2][y2].setNode(nullptr);
    for (auto it = enodes.begin(); it != enodes.end(); ++it) {
        if ((*it).x == a.x && (*it).y == a.y) {
            enodes.erase(it);
            break; // Optional, if you know there is only one element with the value 3
        }
    } //somtimes causes errors
}

void Grid::updateEnodes() {
    enodes.clear();
    for (int i = 1; i < (grid.size() - 1); i++) {
        for (int j = 1; j < (grid[0].size() - 1); j++) {
            if (grid[i][j].getNode() == nullptr) {
                Coords c(i, j);
                enodes.push_back(c);
            }
        }
    }

}

void Grid::move(int x1, int y1, int x2, int y2) {
    if (grid[x1][y1].getType() == squareType::Terminal || grid[x2][y2].getType() == squareType::Terminal) {
        cout << "Cannot cast move on terminal nodes." << endl;
    }
    else if (grid[x1][y1].getType() == squareType::Node && grid[x2][y2].getType() == squareType::Node) {
        if (x1 > grid.size() || y1 > grid[0].size() || x2 > grid.size() || y2 > grid[0].size()) cout << "Trying to Move out of bounds" << endl;
        else  if (grid[x2][y2].getNode() == nullptr) {
            if (grid[x1][y1].getNode() != nullptr) {
                nodeCoords[grid[x1][y1].getNode()->getName()].x = x2;
                nodeCoords[grid[x1][y1].getNode()->getName()].y = y2;
            }
            grid[x2][y2] = grid[x1][y1];
            grid[x1][y1].setNode(nullptr);
            updateEmpties(x1, y1, x2, y2, grid[x1][y1].getType() == squareType::Terminal);
        }
        else if (grid[x1][y1].getNode() == nullptr) {
            nodeCoords[grid[x2][y2].getNode()->getName()].x = x1;
            nodeCoords[grid[x2][y2].getNode()->getName()].y = y1;
            grid[x1][y1] = grid[x2][y2];
            grid[x2][y2].setNode(nullptr);

            updateEmpties(x2, y2, x1, y1, grid[x1][y1].getType() == squareType::Terminal);
        }
        else cout << "Trying to move to none empty Node square" << endl;
    }
}

void Grid::swap(int x1, int y1, int x2, int y2) {
    if (grid[x1][y1].getType() == squareType::Terminal || grid[x2][y2].getType() == squareType::Terminal) {
        cout << "Cannot cast move on terminal nodes." << endl;
    }
    else if (grid[x1][y1].getType() == squareType::Node && grid[x2][y2].getType() == squareType::Node) {
        if (x1 > grid.size() || y1 > grid[0].size() || x2 > grid.size() || y2 > grid[0].size()) cout << "Trying to Swap out of bounds" << endl;
        else if (grid[x2][y2].getNode() == nullptr || grid[x1][y1].getNode() == nullptr) cout << "Trying to swap with a nullptr" << endl;
        else {
            if (grid[x2][y2].getNode() != nullptr) {
                nodeCoords[grid[x2][y2].getNode()->getName()].x = x1;
                nodeCoords[grid[x2][y2].getNode()->getName()].y = y1;
            }
            if (grid[x1][y1].getNode() != nullptr) {
                nodeCoords[grid[x1][y1].getNode()->getName()].x = x2;
                nodeCoords[grid[x1][y1].getNode()->getName()].y = y2;
            }
            square temp = grid[x1][y1];
            grid[x1][y1] = grid[x2][y2];
            grid[x2][y2] = temp;
            //ug.swap(x1, y1, x2, y2);
        }
    }
}

void Grid::smartMutation(int x1, int y1, vector<Bounds> bo, map<string, Net> nets) {
    int ran = rand() % 2;
    int rand_x = 0;
    int rand_y = 0;
    int ran_i = 0;
    Bounds x;
    bool flag = false;
    if (grid[x1][y1].getNode() == nullptr) {
        flag = true;
    }
    if (flag) {
        mutation(x1, y1);
        return;
    }
    for (Bounds b : bo) {
        for (auto net : grid[x1][y1].getNode()->getNets()) {
            if (b.name == net->name) {
                x = b;
                exit;
            }
        }


    }


    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> xdistribution(x.x1, x.x2);
    std::uniform_int_distribution<int> ydistribution(x.y1, x.y2);
    if (x.x2 > 125 || x.y2 > 125) {
        cout << "Crazy Bounds" << endl;
    }
    if (ran % 2 == 0) { //Even will use the swap function
        rand_x = xdistribution(gen);
        rand_y = ydistribution(gen);

        swap(x1, y1, rand_x, rand_y);
    }
    else {  //Odd will use the move function
        for (auto c : enodes) {
            if (c.x < x.x2 && c.x > x.x1 && c.y < x.y2 && c.y > x.y1) {
                if (c.x > 125 || c.y > 125) {
                    cout << "Crazy Coordinate" << endl;
                }
                move(x1, y1, c.x, c.y);
                return;
            }
        }
        ran_i = rand() % enodes.size();
        rand_x = enodes.at(ran_i).x;
        rand_y = enodes.at(ran_i).y;
        move(x1, y1, rand_x, rand_y);
    }
}

void Grid::mutation(int x1, int y1) {
    int ran = rand() % 2;
    int rand_x = 0;
    int rand_y = 0;
    int ran_i = 0;

    if (ran % 2 == 0) { //Even will use the swap function
        rand_x = rand() % grid.size();
        rand_y = rand() % grid[0].size();
        swap(x1, y1, rand_x, rand_y);
    }
    else {  //Odd will use the move function
        ran_i = rand() % enodes.size();
        rand_x = enodes.at(ran_i).x;
        rand_y = enodes.at(ran_i).y;
        move(x1, y1, rand_x, rand_y);
    }
}

square Grid::getSquare(int x, int y) {
    return grid[x][y];
}

void Grid::initialPlacement(const std::map<std::string, Node>& nodes) {
    // Adjust the coordinate system to start from 0,0 if minimum is -33.
    int coordinateShift = 33; // Assuming -33 is the minimum coordinate.
    int maxX = 0, maxY = 0, minX = 0, minY = 0;
    // Containers for edge terminals.
    std::vector<const Node*> topEdge, bottomEdge, leftEdge, rightEdge, isTerminal;
    // For random placement
    std::mt19937 rng{ std::random_device{}() };
    std::set<std::pair<int, int>> occupiedPositions;

    for (auto pair : nodes) {
        if (pair.second.isTerminal()) {
            isTerminal.push_back(&nodes.at(pair.first));
        }
    }

    auto maxElementX = max_element(isTerminal.begin(), isTerminal.end(), [](const Node* a, const Node* b) {//ChatGPT "how to find a max struct in a vector by its int value"
        return a->getX() < b->getX();
        });
    maxX = (*maxElementX)->getX();
    auto maxElementY = max_element(isTerminal.begin(), isTerminal.end(), [](const Node* a, const Node* b) {
        return a->getY() < b->getY();
        });
    maxY = (*maxElementY)->getY();
    auto minElementX = min_element(isTerminal.begin(), isTerminal.end(), [](const Node* a, const Node* b) {
        return a->getX() < b->getX();
        });
    minX = (*minElementX)->getX();
    auto minElementY = min_element(isTerminal.begin(), isTerminal.end(), [](const Node* a, const Node* b) {
        return a->getY() < b->getY();
        });
    minY = (*minElementY)->getY();

    for (auto pair : nodes) {
        Node& node = pair.second;
        if (node.isTerminal()) { // Corrected to use function call syntax
            // Determine the edge for each terminal using getters
            if (node.getY() == maxY) topEdge.push_back(&nodes.at(pair.first)); // Top edge
            else if (node.getY() == minY) bottomEdge.push_back(&nodes.at(pair.first)); // Bottom edge
            else if (node.getX() == minX) leftEdge.push_back(&nodes.at(pair.first)); // Left edge
            else if (node.getX() == maxX) rightEdge.push_back(&nodes.at(pair.first)); // Right edge
        }
    }

    auto distributeTerminals = [&](const std::vector<const Node*>& edgeTerminals, char edge) {
        int numTerminals = edgeTerminals.size();
        int spacing = (edge == 't' || edge == 'b') ? grid[0].size() / (numTerminals + 1)
            : grid.size() / (numTerminals + 1);
        for (int i = 0; i < numTerminals; ++i) {
            int pos = 0;
            if (spacing == 1)  pos = ((i + 1) * spacing);
            else  pos = ((i + 1) * spacing) - 1; //want to start at index 1 no terminals at corners
            int x = 0, y = 0;
            if (edge == 't') { x = pos; y = 0; }
            else if (edge == 'b') { x = pos; y = grid.size() - 1; }
            else if (edge == 'l') { x = 0; y = pos; }
            else if (edge == 'r') { x = grid[0].size() - 1; y = pos; }
            write(y, x, square(squareType::Terminal, edgeTerminals[i]));
        }
    };

    distributeTerminals(topEdge, 't');
    distributeTerminals(bottomEdge, 'b');
    distributeTerminals(leftEdge, 'l');
    distributeTerminals(rightEdge, 'r');

    // Place non-terminal nodes randomly
    std::uniform_int_distribution<int> distX(1, grid.size() - 2), distY(1, grid[0].size() - 2);

    for (auto pair : nodes) {
        auto& node = pair.second;
        if (!node.isTerminal()) {
            bool placed = false;
            while (!placed) {
                int randomX = distX(rng);
                int randomY = distY(rng);
                if (occupiedPositions.find({ randomX, randomY }) == occupiedPositions.end()) {
                    // If position is not occupied, place the node
                    write(randomX, randomY, square(squareType::Node, &nodes.at(node.getName())));
                    occupiedPositions.insert({ randomX, randomY });
                    placed = true;
                }
            }
        }
    }

    for (int i = 1; i < grid.size() - 1; i++) {
        for (int j = 1; j < grid[0].size() - 1; j++) {
            if (grid[i][j].getNode() == nullptr) {
                Coords c(i, j);
                enodes.push_back(c);
            }
        }
    }
}

float Grid::calcCost(float const w1, float const w2, float const w3, const map<string, Net>& nets, bool& routable, int wireConstraint, vector<Bounds>& bounded) const {
    float totalCost = 0, totalLength = 0, overlapCount = 0, critCost = 0;
    bounded.clear(); // Clear previous bounds
    cout << "Calculating cost with weights w1 = " << w1 << ", w2 = " << w2 << ",w3 = " << w3 << endl;
    int gridArea = grid.size() * grid[0].size(); // Calculate grid area once

    // Process each net to calculate wirelength and critical net cost
    for (const auto& [netName, net] : nets) {
        if (net.Nodes.empty()) continue;

        int xmin = INT_MAX, xmax = INT_MIN, ymin = INT_MAX, ymax = INT_MIN;

        // Determine the bounding box for each net
        for (const auto* nodePtr : net.Nodes) {
            const auto& coords = nodeCoords.at(nodePtr->getName());
            xmin = std::min(xmin, coords.x);
            xmax = std::max(xmax, coords.x);
            ymin = std::min(ymin, coords.y);
            ymax = std::max(ymax, coords.y);
        }

        // Calculate total wirelength for the net
        totalLength += (xmax - xmin) + (ymax - ymin);

        // Add additional cost for critical nets
        if (net.isCritical) {
            critCost += ((xmax - xmin) + (ymax - ymin)) / 2; // 50% additional cost for critical nets
        }

        Bounds bounds{ netName, xmin, ymin, xmax, ymax };
        bounded.push_back(bounds);
    }

    // Calculate overlaps and update routability
    overlapCount = calculateOverlaps(bounded, wireConstraint, routable);

    // Normalize cost components
    float normalizedLength = static_cast<float>(totalLength) / gridArea;
    float normalizedOverlap = overlapCount / (nets.size() * gridArea);
    float normalizedCritCost = critCost / (nets.size() * gridArea);

    // Compute total cost
    totalCost = w1 * normalizedLength + w2 * normalizedOverlap + w3 * normalizedCritCost;

    std::cout << "Total Cost:" << totalCost << ", Routable: " << (routable ? "Yes" : "No") << std::endl;
    return totalCost;
}

float Grid::calculateOverlaps(const vector<Bounds>& bounds, int wireConstraint, bool& routable) const {
    float overlapCount = 0;
    routable = true;

    for (size_t i = 0; i < bounds.size(); ++i) {
        for (size_t j = i + 1; j < bounds.size(); ++j) {
            int x_overlap = std::max(0, std::min(bounds[i].x2, bounds[j].x2) - std::max(bounds[i].x1, bounds[j].x1));
            int y_overlap = std::max(0, std::min(bounds[i].y2, bounds[j].y2) - std::max(bounds[i].y1, bounds[j].y1));
            float overlapArea = x_overlap * y_overlap;

            if (overlapArea > 0) {
                overlapCount += overlapArea;
                if (--wireConstraint < 0) routable = false;
            }
        }
    }

    return overlapCount;
}


float updateCost(float const w1, float const w2, float const w3, bool& routable, int wireConstraint, vector<Bounds>& bounded, bool isSwap, int x1, int x2, int y1, int y2) {
    return 0;
}

void Grid::placeNode(int x, int y, const Node* node) {
    if (x >= 0 && x < grid.size() && y >= 0 && y < grid[0].size()) {
        grid[x][y].setNode(node);
        Coords c(x, y);
        nodeCoords[node->getName()] = c;
    }
    else {
        // Handle the error: position out of bounds
        std::cout << "Error: Position (" << x << ", " << y << ") is out of bounds for placing a node.\n";
    }
}

bool Grid::isNodePlaced(const Node* node) const {
    for (vector<square> row : grid) {
        for (auto& square : row) {
            if (square.getNode() == node) return true;
        }
    }
    return false;
}

int Grid::getGridY() {
    return grid.size();
}

int Grid::getGridX() {
    return grid.at(0).size();
}

vector<vector<square>> Grid::getGrid() {
    return grid;
}

map<string, Coords> Grid::getCoords() {
    return nodeCoords;
}