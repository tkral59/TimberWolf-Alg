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
    if (x >= 0 && x < grid.size() && y >= 0 && y < grid[x].size()) {
        grid[x][y] = s;
        if (2*x < ug.grid.size() && 2*y < ug.grid[0].size()) { // Assuming ug.grid is public; adjust if it's private
            ug.write(x * 2, y * 2, s);
        }
    } else {
        // Handle out-of-bounds access appropriately
    }
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

/*void Grid::initialPlacement(const std::map<std::string, Node>& nodes) {
    // Find the scaling factors based on the grid size and maximum coordinates in the .pl file.
    int maxX = 0, maxY = 0;
    for (const auto& node : nodes) {
        maxX = std::max(maxX, node.second.getX());
        maxY = std::max(maxY, node.second.getY());
    }
    
    double scaleX = static_cast<double>(grid.size() - 1) / maxX;
    double scaleY = static_cast<double>(grid[0].size() - 1) / maxY;

    // Place terminal nodes on the perimeter based on their scaled coordinates.
    for (const auto& node : nodes) {
        if (node.second.isTerminal()) {
            // Scale the node's coordinates
            int gridX = static_cast<int>(std::round(node.second.getX() * scaleX));
            int gridY = static_cast<int>(std::round(node.second.getY() * scaleY));

            // Determine the closest edge for each terminal node
            int edgeX = gridX, edgeY = gridY;
            if (gridX <= gridY && gridX <= grid.size() - 1 - gridX) edgeY = 0; // Left edge
            else if (gridY <= gridX && gridY <= grid[0].size() - 1 - gridY) edgeX = 0; // Top edge
            else if (grid.size() - 1 - gridX < gridY) edgeY = grid[0].size() - 1; // Right edge
            else edgeX = grid.size() - 1; // Bottom edge

            // Place the terminal on the closest edge
            write(edgeX, edgeY, square(squareType::Terminal, &node.second));
        }
    }

    // Sequentially place non-terminal nodes in the remaining spaces
    for (const auto& node : nodes) {
        if (!node.second.isTerminal()) {
            // Find an empty space for the non-terminal node
            for (int i = 0; i < grid.size(); ++i) {
                for (int j = 0; j < grid[0].size(); ++j) {
                    if (grid[i][j].getType() == squareType::Routing) { // Assuming Routing type means the space is empty
                        write(i, j, square(squareType::Node, &node.second));
                        goto nextNode;
                    }
                }
            }
            nextNode:;
        }
    }
    
}
*/

void Grid::initialPlacement(const std::map<std::string, Node>& nodes) {
    // Containers for terminals by edge
    std::vector<Node*> topEdge, bottomEdge, leftEdge, rightEdge;

    // Assuming maxX and maxY have been determined from the .pl file
    int maxX = 0, maxY = 0;
    for (const auto& [name, node] : nodes) {
        maxX = std::max(maxX, node.x);
        maxY = std::max(maxY, node.y);
        if (node.isTerminal) {
            // Determine the edge for each terminal
            if (node.y == maxY) topEdge.push_back(const_cast<Node*>(&node)); // Top edge
            else if (node.y == 0) bottomEdge.push_back(const_cast<Node*>(&node)); // Bottom edge
            else if (node.x == 0) leftEdge.push_back(const_cast<Node*>(&node)); // Left edge
            else if (node.x == maxX) rightEdge.push_back(const_cast<Node*>(&node)); // Right edge
        }
    }

    // Helper function to evenly place terminals on their edge
    auto distributeTerminals = [&](std::vector<Node*>& edgeTerminals, char edge) {
        int numTerminals = edgeTerminals.size();
        int spacing = (edge == 't' || edge == 'b') ? grid[0].size() / (numTerminals + 1) 
                                                    : grid.size() / (numTerminals + 1);

        for (int i = 0; i < numTerminals; ++i) {
            int pos = (i + 1) * spacing;
            int x = 0, y = 0;
            switch (edge) {
                case 't': x = pos; y = 0; break;
                case 'b': x = pos; y = grid.size() - 1; break;
                case 'l': x = 0; y = pos; break;
                case 'r': x = grid[0].size() - 1; y = pos; break;
            }
            write(x, y, square(squareType::Terminal, edgeTerminals[i]));
        }
    };

    // Evenly distribute terminals along each edge
    distributeTerminals(topEdge, 't');
    distributeTerminals(bottomEdge, 'b');
    distributeTerminals(leftEdge, 'l');
    distributeTerminals(rightEdge, 'r');

    // Sequentially place non-terminal nodes in the first available spot
    for (const auto& [name, node] : nodes) {
        if (!node.isTerminal) {
            for (int i = 0; i < grid.size(); ++i) {
                for (int j = 0; j < grid[0].size(); ++j) {
                    if (grid[i][j].getType() == squareType::Routing) { // Assuming Routing means empty
                        write(i, j, square(squareType::Node, const_cast<Node*>(&node)));
                        goto nextNode; // Proceed to the next node after placement
                    }
                }
            }
        }
        nextNode: ; // Label to continue with the next non-terminal node
    }
}



float Grid::calcCost(float const w1, float const w2, map<string, Net> const nets, bool& routable, int wireConstraint, vector<Bounds>& bounded) const {
    float totalCost = 0, totalLength = 0, overlapCount = 0;
    //vector<Bounds> bounded;
    bounded.clear();//incase bounded already populated
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

void Grid::placeNode(int x, int y, const Node* node) {
    if (x >= 0 && x < grid.size() && y >= 0 && y < grid[0].size()) {
        grid[x][y].setNode(node);
    }
    else {
        // Handle the error: position out of bounds
        std::cerr << "Error: Position (" << x << ", " << y << ") is out of bounds for placing a node.\n";
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

int Grid::getGridX() {
    return grid.size();
}

int Grid::getGridY() {
    return grid.at(0).size();
}
int Grid::getGridSize() const {
    return grid.size(); // Assuming 'grid' is a container like std::vector
}
