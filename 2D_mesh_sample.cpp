#include <iostream>
#include <vector>
#include <iomanip>
#include <string>

// A simple struct to define the rectangular domain
struct Domain {
    double left, right, bottom, top;
};

// A struct to hold the type of boundary condition for each side
struct BoundaryTypes {
    int bottom; // e.g., -1 for Dirichlet
    int right;  // e.g., -2 for Neumann
    int top;
    int left;
};

/**
 * @brief Generates a 2D mesh and its boundary information based on column-major indexing.
 * 
 * @param P Output, stores node coordinates [x0, y0, x1, y1, ...]
 * @param T Output, stores element connectivity [[n0, n1, n2], ...]
 * @param boundary_edges Output, stores boundary edge data. Each inner vector is 
 *                       {type, element_idx, first_node_idx, second_node_idx}.
 * @param domain The rectangular domain's boundaries.
 * @param N1 Number of intervals in the x-direction.
 * @param N2 Number of intervals in the y-direction.
 * @param bc_types The types of boundary conditions to assign to each of the four sides.
 */
void generateMeshWithBoundaries(
    std::vector<double>& P, 
    std::vector<std::vector<int>>& T,
    std::vector<std::vector<int>>& boundary_edges,
    const Domain& domain, 
    int N1, int N2,
    const BoundaryTypes& bc_types) 
{
    // --- 1. Clear old data ---
    P.clear();
    T.clear();
    boundary_edges.clear();

    if (N1 < 1 || N2 < 1) {
        std::cerr << "Interval counts N1 and N2 must be at least 1." << std::endl;
        return;
    }

    // --- 2. Generate Node Coordinates (P vector) ---
    const int Nx = N1 + 1;
    const int Ny = N2 + 1;
    P.reserve(Nx * Ny * 2);
    const double hx = (domain.right - domain.left) / N1;
    const double hy = (domain.top - domain.bottom) / N2;

    // Column-Major order: Outer loop for columns (x), inner for rows (y)
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            P.push_back(domain.left + i * hx);
            P.push_back(domain.bottom + j * hy);
        }
    }

    // --- 3. Generate Triangular Elements (T vector) ---
    T.reserve(2 * N1 * N2);
    auto node_2d_to_1d = [&](int i, int j) { return i * Ny + j; };

    // Column-Major order for elements
    for (int ie = 0; ie < N1; ++ie) {
        for (int je = 0; je < N2; ++je) {
            int idx_BL = node_2d_to_1d(ie, je);
            int idx_BR = node_2d_to_1d(ie + 1, je);
            int idx_TL = node_2d_to_1d(ie, je + 1);
            int idx_TR = node_2d_to_1d(ie + 1, je + 1);

            // Lower triangle (直角在左下角): 第一个局部编号在左下角 BL
            T.push_back({idx_BL, idx_BR, idx_TL}); 
            // Upper triangle (直角在右上角): 第一个局部编号在左上角 TL
            T.push_back({idx_TL, idx_BR, idx_TR}); 
        }
    }

    // --- 4. Generate Boundary Edges ---
    boundary_edges.reserve(2 * N1 + 2 * N2);
    auto element_2d_to_1d = [&](int ie, int je, bool is_upper) {
        return ie * (2 * N2) + 2 * je + (is_upper ? 1 : 0);
    };

    // BOTTOM boundary (left to right)
    for (int ie = 0; ie < N1; ++ie) {
        int node1 = node_2d_to_1d(ie, 0);
        int node2 = node_2d_to_1d(ie + 1, 0);
        int elem_idx = element_2d_to_1d(ie, 0, false); // Belongs to the lower triangle
        boundary_edges.push_back({bc_types.bottom, elem_idx, node1, node2});
    }

    // RIGHT boundary (bottom to top)
    for (int je = 0; je < N2; ++je) {
        int node1 = node_2d_to_1d(N1, je);
        int node2 = node_2d_to_1d(N1, je + 1);
        int elem_idx = element_2d_to_1d(N1 - 1, je, true); // Belongs to the upper triangle of the last column of cells
        boundary_edges.push_back({bc_types.right, elem_idx, node1, node2});
    }

    // TOP boundary (right to left)
    for (int ie = N1 - 1; ie >= 0; --ie) {
        int node1 = node_2d_to_1d(ie + 1, N2);
        int node2 = node_2d_to_1d(ie, N2);
        int elem_idx = element_2d_to_1d(ie, N2 - 1, true); // Belongs to the upper triangle of the top row of cells
        boundary_edges.push_back({bc_types.top, elem_idx, node1, node2});
    }

    // LEFT boundary (top to bottom)
    for (int je = N2 - 1; je >= 0; --je) {
        int node1 = node_2d_to_1d(0, je + 1);
        int node2 = node_2d_to_1d(0, je);
        int elem_idx = element_2d_to_1d(0, je, false); // Belongs to the lower triangle of the first column of cells
        boundary_edges.push_back({bc_types.left, elem_idx, node1, node2});
    }
}


// --- Main function to demonstrate and verify ---
int main() {
    std::vector<double> P;
    std::vector<std::vector<int>> T;
    std::vector<std::vector<int>> boundary_edges;

    // Recreate the slide's example: N1 = 2, N2 = 2
    int N1 = 2, N2 = 2;
    Domain domain = {0.0, 1.0, 0.0, 1.0}; // A standard unit square

    // Define boundary types: e.g., Dirichlet on top/bottom, Neumann on sides
    BoundaryTypes bc_types = {-1, -2, -1, -2};

    generateMeshWithBoundaries(P, T, boundary_edges, domain, N1, N2, bc_types);

    // --- Print Node Info ---
    std::cout << "--- Node Coordinates (P) ---" << std::endl;
    std::cout << "Total Nodes: " << P.size() / 2 << std::endl;
    for (size_t k = 0; k < P.size() / 2; ++k) {
        std::cout << "Node " << std::setw(2) << k << ": (" 
                  << std::fixed << std::setprecision(2) << P[k * 2] << ", " 
                  << std::fixed << std::setprecision(2) << P[k * 2 + 1] << ")" << std::endl;
    }

    // --- Print Element Info ---
    std::cout << "\n--- Element Connectivity (T) ---" << std::endl;
    std::cout << "Total Elements: " << T.size() << std::endl;
    for (size_t e = 0; e < T.size(); ++e) {
        std::cout << "Element " << std::setw(2) << e << ": [" 
                  << T[e][0] << ", " << T[e][1] << ", " << T[e][2] << "]" << std::endl;
    }

    // --- Print Boundary Edge Info ---
    std::cout << "\n--- Boundary Edges ---" << std::endl;
    std::cout << "Total Boundary Edges: " << boundary_edges.size() << std::endl;
    std::cout << std::left << std::setw(12) << "Edge Idx"
              << std::left << std::setw(8) << "Type"
              << std::left << std::setw(12) << "Elem Idx"
              << std::left << std::setw(14) << "First Node"
              << std::left << std::setw(15) << "Second Node" << std::endl;
    std::cout << std::string(60, '-') << std::endl;

    for (size_t k = 0; k < boundary_edges.size(); ++k) {
        std::cout << std::left << std::setw(12) << k 
                  << std::left << std::setw(8) << boundary_edges[k][0]
                  << std::left << std::setw(12) << boundary_edges[k][1]
                  << std::left << std::setw(14) << boundary_edges[k][2]
                  << std::left << std::setw(15) << boundary_edges[k][3] << std::endl;
    }

    return 0;
}