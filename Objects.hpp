#pragma once
#include <iostream>
#include <vector>
#include <map>
#include <set>
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
    Net(const std::string& name) : name(name) {}
};

struct Coords {
    int x, y;
    Coords(int x, int y) : x(x), y(y) {}
};

struct Bounds {
    int x1, x2, y1, y2;
    const Net* net;
};

class Node {
private:
    std::string name;
    std::vector<Net*> nets; //required forward declaration
    int xcoord, ycoord;

public:
    Node();
    Node(std::string name, std::vector<Net*> nets, int xcoord = 0, int ycoord = 0);
    ~Node();
    int getX() const { return xcoord; }
    int getY() const { return ycoord; }
    std::vector<Net*> getNets() const;
    void addNet(Net* net);
    void removeNet(Net* net);
    std::string getName() const;
    void setCoords(int x, int y);
    const bool isTerminal();
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
    //Node* pins[2];
public:
    square();
    square(squareType type, const Node* n = nullptr, int wires = 0);
    void setType(squareType st);
    squareType getType();
    void incWires(); //increment wire count if wire type
    void decWires(); //dec if wire type
    void setNode(Node* n); //set node if terminal or gate
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
    vector<Coords> enodes; //bounds of empty nodes in grid
    vector<Coords> eterms;
public:
    //Grid();
    Grid(const std::map<std::string, Node>& nodes); // Updated constructor
    void write(int x, int y, square s);
    void swap(int x1, int y1, int x2, int y2);
    void move(int x1, int y1, int x2, int y2);
    void mutation(int x1, int y1);
    void initialPlacement(const std::map<std::string, Node>& nodes);
    void placeTerminals(const std::vector<Node*>& terminals);
    void placeNonTerminals(const std::vector<Node*>& nonTerminals);
    square getSquare(int x, int y); //get square with coordinates

    int calcCost(float const w1, float const w2, map<string, Net> const nets, bool& routable, int wireConstraint) const;
    friend class square;
    friend class utilGrid;
};

#endif // DATASTRUCTURES_HPP
