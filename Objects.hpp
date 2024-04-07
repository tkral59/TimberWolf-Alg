#pragma once
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include<memory>
#ifndef DATASTRUCTURES_HPP
#define DATASTRUCTURES_HPP
using namespace std;
// Forward declaration of Node to resolve circular dependency
class Node;

// Full definition of Net
struct Net {
    Net() = default; // Enable default construction
    std::string name; // Name of the net
    std::vector<Node*> Nodes; // Pointers to Nodes connected to this net
    int nodesSize;
    bool isCritical;
    Net(const std::string& name, int w = 0) : name(name), nodesSize(w) {}

};

struct Coords {
    int x, y;
    Coords() = default;
    Coords(int x, int y) : x(x), y(y) {}
    bool operator==(const Coords& other) const {//ChatGPT
        return x == other.x && y == other.y;
    }
    Node* n;
};

struct Bounds {
    int x1, x2, y1, y2;
    const Net* net;
};

class Node {
private:
    std::string name;
    std::vector<Net*> nets;
    int x = 0, y = 0;
    bool isTerminalFlag = false; // Renamed to avoid conflict with isTerminal() function

public:
    Node(); // Default constructor
    Node(std::string name, std::vector<Net*> nets, int x = 0, int y = 0, bool isTerminal = false);
    ~Node(); // Destructor
    int getX() const;
    int getY() const;
    void setXY(int newX, int newY);
    std::vector<Net*> getNets() const;
    void addNet(Net* net);
    void removeNet(Net* net);
    std::string getName() const;
    bool isTerminal() const; // Fixed return type and made const
    // In Node class in Objects.hpp
    void setTerminal(bool terminal);
    float weight; // Weight attribute for nodes

};


enum class squareType {
    Terminal, //pin e.g. p123 ***in this case the square can hold 2 pins e.g. [p1, p2]
    Node,//gate e.g. a123
    Routing // space for wires
};

class square {
private:
    squareType type;
    const Node* node;
    int wires; //count of wires
public:
    square();
    square(squareType type, const Node* n = nullptr, int wires = 0);
    void setType(squareType st);
    squareType getType();
    void incWires(); //increment wire count if wire type
    void decWires(); //dec if wire type
    void setNode(const Node* n); //set node if terminal or gate
    const Node* getNode();
    bool isEmpty();

    friend class Node;
};

class utilGrid {
private:
    vector<vector<square>> grid;
public:
    utilGrid();
    utilGrid(vector<vector<square>> const ogrid);
    void write(int x, int y, square s); //writing a node into a square
    void swap(int x1, int y1, int x2, int y2); //swapping two nodes
    void move(int x1o, int y1o, int x2o, int y2o);//moving a node to an empty space
    friend class square;
    friend class Grid;
};

class Grid {
private:
    vector<vector<square>> grid;
    utilGrid ug;
    vector<Coords> enodes; //coords of empty nodes in grid
    map<string, Coords> nodeCoords;//map name to coordinates

public:
    Grid();
    Grid(const std::map<std::string, Node>& nodes); // Updated constructor
    void write(int x, int y, square s);
    void swap(int x1, int y1, int x2, int y2);
    void move(int x1, int y1, int x2, int y2);
    void mutation(int x1, int y1);
    void smartMutation(int x1, int y1, vector<Bounds> b);
    void initialPlacement(const std::map<std::string, Node>& nodes);
    square getSquare(int x, int y); //get square with coordinates
    float calcCost(float const w1, float const w2, float const w3, map<string, Net> const nets, bool& routable, int wireConstraint, vector<Bounds>& bounded) const;
    float updateCost(float const w1, float const w2, float const w3, bool& routable, int wireConstraint, vector<Bounds>& bounded, bool isSwap, int x1, int x2, int y1, int y2);
    int getGridX();
    int getGridY();
    void updateEmpties(int x1, int y1, int x2, int y2, bool isTerminal);
    // New methods for crossover support
    //int getGridSize() const;
    void placeNode(int x, int y, const Node* node);
    bool isNodePlaced(const Node* node) const;
    int getGridSize() const { return grid.size(); }
    vector<vector<square>> getGrid();
    map<string, Coords> getCoords();

    // If you decide to make crossover a member function
    //static Grid* crossover(const Grid& parent1, const Grid& parent2, const std::map<std::string, Net>& nets);
    friend class square;
    friend class utilGrid;
};
struct Result {
    Grid g;
    float cost;
    bool routable = false;
    vector<Bounds> bounds;
    Result(Grid g, float cost, bool routable, vector<Bounds> bounds) : g(g), cost(cost), routable(routable), bounds(bounds) {}
    Result() : cost(0), routable(false) {}
};

#endif // DATASTRUCTURES_HPP
