#include <iostream>
#include <array>

#include <fstream>
#include <vector>
#include <unordered_map>
#include <sstream>
#include <string>
#include <algorithm>
#include <cmath>
#include <omp.h>
#include <stdexcept>
#include <iomanip>

using namespace std;
using Tetrahedron = std::array<int, 4>;


// Define a 3D point
struct Point {
    double x, y, z;

    // Subtract operator
    Point operator-(const Point& other) const {
        return {x - other.x, y - other.y, z - other.z};
    }

    // Dot product
    double dot(const Point& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    // Cross product
    Point cross(const Point& other) const {
        return {
            y * other.z - z * other.y,
            z * other.x - x * other.z,
            x * other.y - y * other.x
        };
    }
};


std::pair<std::unordered_map<int, Point>, std::vector<Tetrahedron>> parseNastranShortField(const std::string& filePath, Point& minCoords, Point& maxCoords) {
    std::unordered_map<int, Point> nodes;
    std::vector<Tetrahedron> tetrahedra;

    std::ifstream file(filePath);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file");
    }

    minCoords = {std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
    maxCoords = {std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest()};

    std::string line;
    while (std::getline(file, line)) {
        if (line.substr(0, 4) == "GRID") {
            int nodeId = std::stoi(line.substr(8, 8));
            double x = std::stod(line.substr(24, 8));
            double y = std::stod(line.substr(32, 8));
            double z = std::stod(line.substr(40, 8));
            cout <<"x y z "<<x <<", "<<y <<", "<<z<<endl;
            nodes[nodeId] = {x, y, z};

            minCoords.x = std::min(minCoords.x, x);
            minCoords.y = std::min(minCoords.y, y);
            minCoords.z = std::min(minCoords.z, z);

            maxCoords.x = std::max(maxCoords.x, x);
            maxCoords.y = std::max(maxCoords.y, y);
            maxCoords.z = std::max(maxCoords.z, z);
        } else if (line.substr(0, 6) == "CTETRA") {
            Tetrahedron tet;
            tet[0] = std::stoi(line.substr(24, 8));
            tet[1] = std::stoi(line.substr(32, 8));
            tet[2] = std::stoi(line.substr(40, 8));
            tet[3] = std::stoi(line.substr(48, 8));
            tetrahedra.push_back(tet);
            cout << "tetra "<<tet[0]<<", "<<tet[1]<<", "<<tet[2]<<", "<<tet[3]<<endl;
        }
    }

    return {nodes, tetrahedra};
}

/*
////IS ELEMENT IS ORIENTED AND NORMAL ARE ALWAYS OUTER OR INNER
//CHANGE TO TETRAHEDRON
// Compute if a point is inside a tetrahedron
bool isPointInsideTetra(const Point& p, const Point& a, const Point& b, const Point& c, const Point& d) {
    // Compute normals for each face
    Point n1 = (b - a).cross(c - a); // Normal to face ABC
    Point n2 = (c - a).cross(d - a); // Normal to face ACD
    Point n3 = (d - a).cross(b - a); // Normal to face ADB
    Point n4 = (c - b).cross(d - b); // Normal to face BCD

    
    // Check the dot product with the point and a vertex of each face
    bool inside = true;
    inside &= (n1.dot(p - a) <= 0); // Check face ABC
    inside &= (n2.dot(p - a) <= 0); // Check face ACD
    inside &= (n3.dot(p - a) <= 0); // Check face ADB
    inside &= (n4.dot(p - b) <= 0); // Check face BCD

    return inside;
}
*/

bool sameSide(const Point& p, const Point& a, const Point& b, const Point& c, const Point& d) {
    bool ret = false;
    Point n = (b - a).cross(c - a); // Normal to face ABC  
    double dotd = n.dot(d-a);
    double dotp = n.dot(p-a);
    if (dotd>0&&dotp>0 || dotd<0 && dotp<0)
      ret = true;
    
    return ret;
}
//CHANGE TO TETRAHEDRON
// Compute if a point is inside a tetrahedron
bool isPointInsideTetra(const Point& p, const Point& a, const Point& b, const Point& c, const Point& d) {
    // Compute normals for each face
    return sameSide(p,a,b,c,d) &&
           sameSide(p,b,c,d,a) &&
           sameSide(p,c,d,a,b) &&
           sameSide(p,d,a,b,c);
    
}


bool isPointInsideMesh(const Point& point, const std::unordered_map<int, Point>& nodes, const std::vector<Tetrahedron>& tetrahedra) {
    int i=0;
    for (const auto& tet : tetrahedra) {
        //cout << "TETRA "<<i<<endl;
        //std::array<Point, 4> vertices = {nodes.at(tet[0]), nodes.at(tet[1]), nodes.at(tet[2]), nodes.at(tet[3])};
        if (isPointInsideTetra(point, /*vertices*/nodes.at(tet[0]), nodes.at(tet[1]), nodes.at(tet[2]), nodes.at(tet[3]))) {
            return true;
        }
        i++;
    }
    return false;
}

std::vector<Point> generateCartesianGrid(const Point& minCoords, const Point& maxCoords, double spacing) {
    std::vector<Point> grid;
    for (double x = minCoords.x; x <= maxCoords.x; x += spacing) {
        for (double y = minCoords.y; y <= maxCoords.y; y += spacing) {
            for (double z = minCoords.z; z <= maxCoords.z; z += spacing) {
                grid.push_back({x, y, z});
            }
        }
    }
    std::cout << "Generating "<<grid.size() << " points "<<std::endl;
    return grid;
}
/*
void writeSphToKFile(const std::vector<Point>& points, const std::string& outputFile) {
    std::ofstream file(outputFile);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open output file");
    }

    file << "*KEYWORD\n";
    file << "*NODE\n";
    for (size_t i = 0; i < points.size(); ++i) {
        file << i + 1 << "," << points[i].x << "," << points[i].y << "," << points[i].z << "\n";
    }

    file << "*ELEMENT_SPH\n";
    for (size_t i = 0; i < points.size(); ++i) {
        file << i + 1 << "," << i + 1 << "\n";
    }

    file << "*END\n";
}
*/
void writeSphToKFile(const std::vector<Point>& points, const std::string& outputFile) {
    std::ofstream file(outputFile);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open output file");
    }

    file << "*KEYWORD\n";
    file << "*NODE\n";

    for (size_t i = 0; i < points.size(); ++i) {
        file 
            << std::setw(8) << i + 1 // Node ID
            << std::setw(16) << std::fixed << std::setprecision(6) << points[i].x // X-coordinate
            << std::setw(16) << points[i].y // Y-coordinate
            << std::setw(16) << points[i].z // Z-coordinate
            << "\n";
    }

    file << "*ELEMENT_SPH\n";

    for (size_t i = 0; i < points.size(); ++i) {
        file 
            << std::setw(8) << i + 1 // Element ID
            << std::setw(8) << i + 1 // Node ID
            << "\n";
    }

    file << "*END\n";
}

/*
int main() {
    // Tetrahedron vertices
    Point a = {0, 0, 0};
    Point b = {1, 0, 0};
    Point c = {0, 1, 0};
    Point d = {0, 0, 1};

    // Point to check
    Point p = {0.25, 0.25, 0.25}; // Inside
    Point q = {1, 1, 1};          // Outside

    std::cout << "Point P is " << (isPointInsideTetra(p, a, b, c, d) ? "inside" : "outside") << " the tetrahedron.\n";
    std::cout << "Point Q is " << (isPointInsideTetra(q, a, b, c, d) ? "inside" : "outside") << " the tetrahedron.\n";

    return 0;
}*/


int main(int argc, char *argv[]) {

    double spacing = 1.0;
    if (--argc>0){
      spacing = atof(argv[argc]);
    }
    
    cout << "Spacing is "<<spacing<<endl;
    std::string nastranFile = "mesh.bdf";
    Point minCoords, maxCoords;

    auto [nodes, tetrahedra] = parseNastranShortField(nastranFile, minCoords, maxCoords);
    cout << "Nodes: "<<nodes.size()<<", tetras "<<tetrahedra.size()<<endl;
    std::cout << "Bounding box: Min(" << minCoords.x << ", " << minCoords.y << ", " << minCoords.z << "), Max(" << maxCoords.x << ", " << maxCoords.y << ", " << maxCoords.z << ")\n";


    auto grid = generateCartesianGrid(minCoords, maxCoords, spacing);
    cout << "Done "<<endl;
    
    std::vector<Point> insidePoints;
    int ins = 0;
    #pragma omp parallel for
    for (size_t i = 0; i < grid.size(); ++i) {
        if (i%100 == 0)
          cout << "POINT "<<i<<" of "<<grid.size()<<endl;
        //cout << "POINT "<<i<<endl;
        //print(grid[i]);
        if (isPointInsideMesh(grid[i], nodes, tetrahedra)) {
            //cout << "Grid "<<i<<grid[i].x<<", "<<grid[i].y<<endl;
            #pragma omp critical
            insidePoints.push_back(grid[i]);
            ins ++;
        }
        //cout << "Chekcing point "<<i <<", Inside points "<<ins<<endl; 
    }

    std::string outputFile = "sph_particles.k";
    writeSphToKFile(insidePoints, outputFile);

    std::cout << "SPH particles written to " << outputFile << "\n";
    return 0;
}

