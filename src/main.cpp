#include <iostream>
#include <vector>
#include <cstring>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers> 
#include "SymmetricSparseMatrices.h"
#include <functional>

#define REAL double
#define VOID void
extern "C" {
    #include "triangle.h"
}

struct Mesh {
    Eigen::MatrixX2d V;
    Eigen::MatrixX3i F;
};

double get_boundary_p(double x, double y, double t) {
    return (1.0 - std::exp(-2.0 * t)) + 0.5 * std::sin(2 * t);
}

using BoundaryFunc = std::function<double(double x, double y, double t)>;




Mesh generate_mesh() {
    struct triangulateio in, out;
    std::memset(&in, 0, sizeof(in));
    std::memset(&out, 0, sizeof(out));

    int num_boundary_points = 40; 
    double a = 2.0; 
    double b = 1.0; 

    in.numberofpoints = num_boundary_points;
    in.pointlist = (REAL*)malloc(in.numberofpoints * 2 * sizeof(REAL));
    
    in.numberofsegments = num_boundary_points;
    in.segmentlist = (int*)malloc(in.numberofsegments * 2 * sizeof(int));

    for (int i = 0; i < num_boundary_points; ++i) {
        double angle = 2.0 * M_PI * i / num_boundary_points;
        in.pointlist[2 * i] = a * std::cos(angle);     // X
        in.pointlist[2 * i + 1] = b * std::sin(angle); // Y

        in.segmentlist[2 * i] = i;
        in.segmentlist[2 * i + 1] = (i + 1) % num_boundary_points;
    }

    triangulate((char*)"pq30a0.01zQ", &in, &out, nullptr);

    Mesh mesh;
    mesh.V.resize(out.numberofpoints, 2);
    for (int i = 0; i < out.numberofpoints; ++i) {
        mesh.V(i, 0) = out.pointlist[i * 2];
        mesh.V(i, 1) = out.pointlist[i * 2 + 1];
    }
    mesh.F.resize(out.numberoftriangles, 3);
    for (int i = 0; i < out.numberoftriangles; ++i) {
        mesh.F(i, 0) = out.trianglelist[i * 3];
        mesh.F(i, 1) = out.trianglelist[i * 3 + 1];
        mesh.F(i, 2) = out.trianglelist[i * 3 + 2];
    }

    free(in.pointlist); free(in.segmentlist);
    return mesh;
}

void assemble_flux(const Mesh& mesh, double K, SymmetricSparseMatrices& sm) {
    std::fill(sm.jacobian.begin(), sm.jacobian.end(), 0.0);
    std::fill(sm.residual.begin(), sm.residual.end(), 0.0);
    for (int t = 0; t < mesh.F.rows(); ++t) {
        int ids[3] = {mesh.F(t, 0), mesh.F(t, 1), mesh.F(t, 2)};
        Eigen::Vector2d coords[3];
        for(int i=0; i<3; ++i) coords[i] = mesh.V.row(ids[i]);
        
        double area = 0.5 * std::abs(coords[0].x()*(coords[1].y() - coords[2].y()) + coords[1].x()*(coords[2].y() - coords[0].y()) + coords[2].x()*(coords[0].y() - coords[1].y()));
        
        for (int i = 0; i < 3; ++i) {
            int id_i = ids[i], id_j = ids[(i+1)%3];
            Eigen::Vector2d vi = coords[i] - coords[(i+2)%3];
            Eigen::Vector2d vj = coords[(i+1)%3] - coords[(i+2)%3];
            double T_ij = (K * (vi.dot(vj))) / (4.0 * area);
            
            sm.updateValue(id_i, id_j, -T_ij, -T_ij);
            sm.updateValue(id_j, id_i, -T_ij, -T_ij);
            sm.updateValue(id_i, id_i, T_ij, T_ij);
            sm.updateValue(id_j, id_j, T_ij, T_ij);
        }
    }
}

double get_M(double p) { return 0.01 * (p + 0.01 * p * p); }
double get_dM(double p) { return 0.01 * (1.0 + 0.02 * p); }

int main() {
    Mesh mesh = generate_mesh();
    int n = mesh.V.rows();
    std::cout << "Mesh generated with " << n << " nodes." << std::endl;

    std::vector<std::pair<int, int>> edges;
    for (int i = 0; i < mesh.F.rows(); ++i) {
        for(int j=0; j<3; ++j) edges.push_back({mesh.F(i, j), mesh.F(i, (j+1)%3)});
    }
    
    SymmetricSparseMatrices flux_matrices;
    flux_matrices.initialize(n, edges);

    Eigen::VectorXd p = Eigen::VectorXd::Constant(n, 1.0); 
    Eigen::VectorXd p_old = p;
    double dt = 0.1; 
    double total_steps = 200;

    std::vector<bool> is_dirichlet(n, false);
    for (int i = 0; i < n; ++i) {
        if (mesh.V(i, 0) < 1e-9) is_dirichlet[i] = true;
    }
    double well_x = 0.3, well_y = 0.6;
    double well_q = -1.0; 

    for (int step = 0; step < total_steps; ++step) {
        double current_time = step * dt;
        std::cout << "\n--- Time Step " << step + 1 << " (t=" << current_time << ") ---" << std::endl;

        for (int iter = 0; iter < 15; ++iter) {
            assemble_flux(mesh, 1.0, flux_matrices);
            
            Eigen::SparseMatrix<double> A(n, n);
            std::vector<Eigen::Triplet<double>> triplets;
            Eigen::VectorXd b = Eigen::VectorXd::Zero(n);

            for (int i = 0; i < n; ++i) {
                if (is_dirichlet[i]) {
                    triplets.push_back({i, i, 1.0});
                    b[i] = get_boundary_p(mesh.V(i,0), mesh.V(i,1), current_time) - p[i]; 
                }
                else {
                    for (int k = flux_matrices.outerIndexPtr[i]; k < flux_matrices.outerIndexPtr[i+1]; ++k) {
                        int col = flux_matrices.innerIndexPtr[k];
                        triplets.push_back({i, col, flux_matrices.jacobian[k]});
                        b[i] -= flux_matrices.residual[k] * p[col];
                }
                    

                    double acc_der = get_dM(p[i]) / dt;
                    double acc_res = (get_M(p[i]) - get_M(p_old[i])) / dt;
                    triplets.push_back({i, i, acc_der});
                    b[i] -= acc_res;

                    double dist = std::sqrt(std::pow(mesh.V(i,0)-well_x, 2) + std::pow(mesh.V(i,1)-well_y, 2));
                    if (dist < 0.1) {
                        b[i] += well_q; 
                    }
                }
            }

            A.setFromTriplets(triplets.begin(), triplets.end());
            Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
            solver.compute(A);
            
            if(solver.info() != Eigen::Success) break;
            
            Eigen::VectorXd delta_p = solver.solve(b);
            p += delta_p;
            if (b.norm() < 1e-7) break;
        }
        p_old = p;
    }

   
    std::ofstream outFile("solution.csv");
    outFile << "x,y,p\n";
    for (int i = 0; i < n; ++i) {
        outFile << mesh.V(i, 0) << "," << mesh.V(i, 1) << "," << p[i] << "\n";
    } 
    outFile.close();
    
    std::cout << "\nDone! Saved to pressure_data.csv. Ready for Python plotting!" << std::endl;
    return 0;
}