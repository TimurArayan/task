#pragma once
#include <vector>
#include <Eigen/Sparse>

class SymmetricSparseMatrices {
public:
    std::vector<int> outerIndexPtr;
    std::vector<int> innerIndexPtr;

    std::vector<double> jacobian;
    std::vector<double> residual;

    int rows = 0;

   
    void initialize(int numNodes, const std::vector<std::pair<int, int>>& edges);
    
    void updateValue(int row, int col, double jac_val, double res_val);
};