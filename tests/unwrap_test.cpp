#include <iostream>
#include <vector>

#include "unwrap.h"
#include "vtkWriter.h"

typedef std::array<int, 3> triangle;

int main(int argc, char **argv) {
    // cerate vertices
    paraFEM::Vector3 p0(-1, -1, 0);
    paraFEM::Vector3 p1(1, -1, 0);
    paraFEM::Vector3 p2(1, 1, 0);
    paraFEM::Vector3 p3(-1, 1, 0);
    paraFEM::Vector3 p4(0, 0, 1);

    std::vector<paraFEM::Vector3> vertices = {p0, p1, p2, p3, p4};

    // create triangles
    std::vector<triangle> triangles;
    triangles.push_back(triangle{ {0, 1, 4} });
    triangles.push_back(triangle{ {1, 2, 4} });
    triangles.push_back(triangle{ {2, 3, 4} });
    triangles.push_back(triangle{ {3, 0, 4} });

    paraFEM::VtkWriter writer = paraFEM::VtkWriter("/tmp/paraFEM/unwrap");

    paraFEM::LscmRelax flattener(vertices, triangles, std::vector<int>());
    for (auto point: flattener.get_flat_vertices())
        std::cout << point << std::endl;
    flattener.lscm();
    for (auto point: flattener.get_flat_vertices())
        std::cout << point << std::endl;

}
