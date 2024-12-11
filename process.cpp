#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <unordered_map>
#include <sstream>
#include <string>
#include <algorithm>
#include <cmath>
#include <omp.h>

struct Point {
    double x, y, z;
};

using Tetrahedron = std::array<int, 4>;
using namespace std;

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
double signedVolume(const Point& a, const Point& b, const Point& c, const Point& d) {
    return ((a.x * (b.y * c.z - c.y * b.z) - a.y * (b.x * c.z - c.x * b.z) + a.z * (b.x * c.y - c.x * b.y)) -
            (a.x * (b.y * d.z - d.y * b.z) - a.y * (b.x * d.z - d.x * b.z) + a.z * (b.x * d.y - d.x * b.y))) / 6.0;
}
*/


// Compute cross product of two 3D vectors
Point crossProduct(const Point& u, const Point& v) {
    return {
        u.y * v.z - u.z * v.y,
        u.z * v.x - u.x * v.z,
        u.x * v.y - u.y * v.x
    };
}

// Compute dot product of two 3D vectors
double dotProduct(const Point& u, const Point& v) {
    return u.x * v.x + u.y * v.y + u.z * v.z;
}

/*
// Compute signed volume of a tetrahedron
double signedVolume(const Point& a, const Point& b, const Point& c, const Point& d) {
    Point ab = {b.x - a.x, b.y - a.y, b.z - a.z};
    Point ac = {c.x - a.x, c.y - a.y, c.z - a.z};
    Point ad = {d.x - a.x, d.y - a.y, d.z - a.z};

    Point cross = crossProduct(ac, ad);
    double volume = dotProduct(ab, cross) / 6.0;
    return volume;
}
*/

double signedVolume(const Point& a, const Point& b, const Point& c, const Point& d) {
    Point ab = {b.x - a.x, b.y - a.y, b.z - a.z};
    Point ac = {c.x - a.x, c.y - a.y, c.z - a.z};
    Point ad = {d.x - a.x, d.y - a.y, d.z - a.z};

    double volume = 1.0/6.0 * (ab.x*(ac.y*ad.z-ac.z*ad.y)-ab.y*(ac.x*ad.z-ac.z*ad.x)+ab.z*(ac.x*ad.y-ac.y*ad.x));
    
    //Point cross = crossProduct(ac, ad);
    //double volume = dotProduct(ab, cross) / 6.0;
    
    return volume;
    
}


void print(const Point & p){
  cout << p.x<<", "<<p.y<<", "<<p.z<<endl;
  }
bool isPointInsideTetrahedron(const Point& p, const std::array<Point, 4>& vertices) {
    for (auto &v:vertices )
      print (v);
      
    double vABCD = signedVolume(vertices[0], vertices[1], vertices[2], vertices[3]);
    double vPBCD = signedVolume(p, vertices[1], vertices[2], vertices[3]);
    double vPADC = signedVolume(vertices[0], p, vertices[3], vertices[2]);
    double vPABD = signedVolume(vertices[0], vertices[1], p, vertices[3]);
    double vPABC = signedVolume(vertices[0], vertices[1], vertices[2],p);
    cout << "ABS "<< std::fabs(vPBCD + vPADC + vPABD + vPABC - vABCD)<<endl;
    return std::fabs(vPBCD + vPADC + vPABD + vPABC - vABCD) < 1e-1;
}

bool isPointInsideMesh(const Point& point, const std::unordered_map<int, Point>& nodes, const std::vector<Tetrahedron>& tetrahedra) {
    int i=0;
    for (const auto& tet : tetrahedra) {
        cout << "TETRA "<<i<<endl;
        std::array<Point, 4> vertices = {nodes.at(tet[0]), nodes.at(tet[1]), nodes.at(tet[2]), nodes.at(tet[3])};
        if (isPointInsideTetrahedron(point, vertices)) {
            cout << "FOUND"<<endl;
            return true;
            cout << "FOUND"<<endl;
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
        cout << "POINT "<<i<<endl;
        print(grid[i]);
        if (isPointInsideMesh(grid[i], nodes, tetrahedra)) {
            cout << "Grid "<<i<<grid[i].x<<", "<<grid[i].y<<endl;
            #pragma omp critical
            insidePoints.push_back(grid[i]);
            ins ++;
        }
        cout << "Chekcing point "<<i <<", Inside points "<<ins<<endl; 
    }

    std::string outputFile = "sph_particles.k";
    writeSphToKFile(insidePoints, outputFile);

    std::cout << "SPH particles written to " << outputFile << "\n";
    return 0;
}

