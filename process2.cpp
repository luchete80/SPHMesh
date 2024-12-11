#include <iostream>
#include <array>
#include <cmath>
#include <vector>
using namespace std;
// Define a 3D point
struct Point {
    double x, y, z;
};


/*
// Calculate the signed volume of a tetrahedron formed by four points
double signedVolume(const Point& a, const Point& b, const Point& c, const Point& d) {
    return (a.x * (b.y * c.z - c.y * b.z) - a.y * (b.x * c.z - c.x * b.z) + a.z * (b.x * c.y - c.x * b.y) -
            d.x * (b.y * c.z - c.y * b.z) + d.y * (b.x * c.z - c.x * b.z) - d.z * (b.x * c.y - c.x * b.y)) / 6.0;
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

// Compute signed volume of a tetrahedron
double signedVolume(const Point& a, const Point& b, const Point& c, const Point& d) {
    Point ab = {b.x - a.x, b.y - a.y, b.z - a.z};
    Point ac = {c.x - a.x, c.y - a.y, c.z - a.z};
    Point ad = {d.x - a.x, d.y - a.y, d.z - a.z};

    double volume = 1.0/6.0 * (ab.x*(ac.y*ad.z-ac.z*ad.y)-ab.y*(ac.x*ad.z-ac.z*ad.x)+ab.z*(ac.x*ad.y-ac.y*ad.x));
    /*
    Point cross = crossProduct(ac, ad);
    double volume = dotProduct(ab, cross) / 6.0;
    */
    return volume;
    
}

// Check if a point is inside a tetrahedron
bool isPointInsideTetrahedron(const Point& p, const std::array<Point, 4>& vertices) {
    double vOriginal = signedVolume(vertices[0], vertices[1], vertices[2], vertices[3]);
    double v1 = signedVolume(p, vertices[1], vertices[2], vertices[3]);
    double v2 = signedVolume(vertices[0], p, vertices[2], vertices[3]);
    double v3 = signedVolume(vertices[0], vertices[1], p, vertices[3]);
    double v4 = signedVolume(vertices[0], vertices[1], vertices[2], p);

    double vSum = v1 + v2 + v3 + v4;
    cout << "voriginal "<<vOriginal<<endl;
    cout << "vsum" <<vSum<<endl;
    // Check if the sum of volumes matches and all volumes have the same sign
    return std::fabs(vSum - vOriginal) < 1e-6 &&
           ((v1 > 0 && v2 > 0 && v3 > 0 && v4 > 0) || (v1 < 0 && v2 < 0 && v3 < 0 && v4 < 0));
}

int main() {
    // Tetrahedron vertices
    std::array<Point, 4> vertices = {
        Point{0, 0, 0},
        Point{0, 1, 0},
        Point{1, 0, 0},
        Point{0, 0, 1}
    };

    // Define the cube range
    double cubeMin = 0.0;
    double cubeMax = 1.0;
    double step = 0.1; // Grid spacing

    std::vector<Point> insidePoints;

    // Iterate through all points in the cube
    for (double x = cubeMin; x <= cubeMax; x += step) {
        for (double y = cubeMin; y <= cubeMax; y += step) {
            for (double z = cubeMin; z <= cubeMax; z += step) {
                Point p = {x, y, z};
                if (isPointInsideTetrahedron(p, vertices)) {
                    insidePoints.push_back(p);
                    std::cout << "FOUND "<<std::endl;
                }
            }
        }
    }

    // Output the points inside the tetrahedron
    std::cout << "Points inside the tetrahedron:\n";
    for (const auto& point : insidePoints) {
        std::cout << "(" << point.x << ", " << point.y << ", " << point.z << ")\n";
    }

    return 0;
}
