#include "SymmetricSparseMatrices.h"
#include <set>
#include <algorithm>
#include <stdexcept>

void SymmetricSparseMatrices::initialize(int numNodes, const std::vector<std::pair<int, int>>& edges) {
    this->rows = numNodes;
    std::vector<std::set<int>> adjacency(numNodes);
    
    for (const auto& edge : edges) {
        int u = edge.first;
        int v = edge.second;
        if (u >= numNodes || v >= numNodes) continue;
        adjacency[u].insert(v);
        adjacency[v].insert(u);
        adjacency[u].insert(u);
        adjacency[v].insert(v);
    }

    outerIndexPtr.clear();
    innerIndexPtr.clear();
    outerIndexPtr.push_back(0);

    for (int i = 0; i < numNodes; ++i) {
        for (int neighbor : adjacency[i]) {
            innerIndexPtr.push_back(neighbor);
        }
        outerIndexPtr.push_back(static_cast<int>(innerIndexPtr.size()));
    }

    int nnz = static_cast<int>(innerIndexPtr.size());
    jacobian.assign(nnz, 0.0);
    residual.assign(nnz, 0.0);
}


void SymmetricSparseMatrices::updateValue(int row, int col, double jac_val, double res_val) {
    for (int k = outerIndexPtr[row]; k < outerIndexPtr[row + 1]; ++k) {
        if (innerIndexPtr[k] == col) {
            jacobian[k] += jac_val; 
            residual[k] += res_val; 
            return;
        }
    }
}