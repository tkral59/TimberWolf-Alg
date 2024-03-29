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

const bool Node::isTerminal() {
    return name[0] == 'p';
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

void Grid::move(int x1, int y1, int x2, int y2) {
    if ((grid[x1][y1].getType() == squareType::Node && grid[x2][y2].getType() == squareType::Node) || (grid[x1][y1].getType() == squareType::Terminal && grid[x2][y2].getType() == squareType::Terminal)) { // if one is node and the other is empty
        if (x1 > grid.size() || y1 > grid[0].size() || x2 > grid.size() || y2 > grid[0].size()) cout << "movement dims outside of grid." << endl;
        else {
            if (grid[x2][y2].getNode() == nullptr) { grid[x2][y2] = grid[x1][y1]; }
            else if (grid[x1][y1].getNode() == nullptr) { grid[x1][y1] = grid[x2][y2]; }

            ug.move(x1, y1, x2, y2);
        }
    }
    else if (grid[x1][y1].getType() == squareType::Node || grid[x2][y2].getType() == squareType::Node || grid[x1][y1].getType() == squareType::Terminal || grid[x2][y2].getType() == squareType::Terminal) cout << "Error: Trying to move non-matching types." << endl;
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
        if (getSquare(x1, y1).getType() == squareType::Terminal) {
            ran_i = rand() % eterms.size();
            rand_x = eterms.at(ran_i).x;
            rand_y = eterms.at(ran_i).y;
            move(x1, y1, rand_x, rand_y);
        }
        else {
            ran_i = rand() % enodes.size();
            rand_x = enodes.at(ran_i).x;
            rand_y = enodes.at(ran_i).y;
            move(x1, y1, rand_x, rand_y);
        }
    }
}

square Grid::getSquare(int x, int y) {
    return grid[x][y];
}

void Grid::initialPlacement(const std::map<std::string, Node>& nodes) {
    std::vector<Node*> terminals, nonTerminals;
    // Separate terminals and non-terminals
    for (auto node : nodes) {
        if (node.second.isTerminal()) terminals.push_back(&node.second);
        else nonTerminals.push_back(&node.second);
    }

    // Shuffle terminals for random perimeter placement
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(terminals.begin(), terminals.end(), g);

    // Place terminals on the perimeter
    int perimeterCount = 0, gridSize = grid.size();
    int maxPerimeterPositions = (gridSize - 1) * 4; // Calculate available perimeter positions
    int terminaliter = 0;
    while (perimeterCount < maxPerimeterPositions) {
        // Calculate position (x,y) based on perimeterCount
        int x = perimeterCount % (gridSize - 1);
        int y = perimeterCount / (gridSize - 1);
        if (y == 1) x = gridSize - 1;
        else if (y == 2) x = gridSize - 1 - x;
        else if (y == 3) y = gridSize - 1;
        else y = 0;

        if (terminaliter < terminals.size()) {
            write(x, y, square(squareType::Terminal, terminals.at(terminaliter))); // Adjust for your implementation
            terminaliter++;
        }
        else {
            write(x, y, square(squareType::Terminal, nullptr));
            Coords ec = Coords(x, y);
            eterms.push_back(ec);
        }
        perimeterCount++;
    }

    // Shuffle and place non-terminals within the grid's interior
    std::shuffle(nonTerminals.begin(), nonTerminals.end(), g);
    int ntiter = 0;
    for (int i = 1; i < gridSize - 1; i++) {
        for (int j = 1; j < gridSize - 1; j++) {
            if (ntiter < nonTerminals.size()) {
                write(i, j, square(squareType::Node, nonTerminals.at(ntiter)));
                ntiter++;
            }
            else { // Assuming Routing indicates an empty spot
                write(i, j, square(squareType::Node, nullptr)); // Adjust for your implementation
                Coords ec = Coords(i, j);
                enodes.push_back(ec);
            }
        }
    }

}

void Grid::placeTerminals(const std::vector<Node*>& terminals) {
    int perimeterIndex = 0;
    // Calculating the perimeter's length excluding the corners twice
    int perimeterLength = 2 * (grid.size() + grid[0].size()) - 4;

    for (Node* terminal : terminals) {
        if (perimeterIndex >= perimeterLength) {
            std::cerr << "Error: More terminals than perimeter slots available.\n";
            break; // Or handle this case as needed
        }

        int x, y;
        // Translate perimeterIndex to (x, y) coordinates on the perimeter
        if (perimeterIndex < grid.size()) { // Top row
            x = perimeterIndex;
            y = 0;
        }
        else if (perimeterIndex < grid.size() + grid[0].size() - 2) { // Right column
            x = grid.size() - 1;
            y = perimeterIndex - grid.size() + 1;
        }
        else if (perimeterIndex < 2 * grid.size() + grid[0].size() - 4) { // Bottom row
            x = 2 * grid.size() + grid[0].size() - 5 - perimeterIndex;
            y = grid[0].size() - 1;
        }
        else { // Left column
            x = 0;
            y = perimeterLength - perimeterIndex;
        }

        grid[y][x] = square(squareType::Terminal, terminal); // Assuming direct assignment is valid
        perimeterIndex++;
    }
}

void Grid::placeNonTerminals(const std::vector<Node*>& nonTerminals) {
    std::random_device rd;
    std::mt19937 g(rd());

    // Gather all internal grid positions
    std::vector<std::pair<int, int>> positions;
    for (int y = 1; y < grid.size() - 1; y++) {
        for (int x = 1; x < grid[0].size() - 1; x++) {
            // Assuming squareType::Routing implies an empty square
            if (grid[y][x].getType() == squareType::Routing) {
                positions.emplace_back(x, y);
            }
        }
    }

    // Shuffle positions for random placement
    std::shuffle(positions.begin(), positions.end(), g);

    for (size_t i = 0; i < nonTerminals.size() && i < positions.size(); i++) {
        int x = positions[i].first, y = positions[i].second;
        grid[y][x] = square(squareType::Node, nonTerminals[i]); // Assuming direct assignment is valid
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
