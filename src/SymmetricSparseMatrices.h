#pragma once
#include <vector>
#include <Eigen/Sparse>

class SymmetricSparseMatrices {
public:
    // 1. Общая структура разреженности (ТЗ: массивы outerIndexPtr, innerIndexPtr)
    std::vector<int> outerIndexPtr;
    std::vector<int> innerIndexPtr;

    // 2. Два отдельных массива значений: якобиан и невязка
    std::vector<double> jacobian;
    std::vector<double> residual;

    int rows = 0;

    // Инициализация структуры разреженности по списку ребер сетки
    void initialize(int numNodes, const std::vector<std::pair<int, int>>& edges);

    // Метод для эффективного обновления значений в обеих матрицах одновременно
    // Мы будем использовать его при сборке потокового члена
    void updateValue(int row, int col, double jac_val, double res_val);
};