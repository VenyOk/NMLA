#include <iostream>
#include <vector>
#include <cmath>
#include <cblas.h>  // Для BLAS функций

class CholeskySolver {
private:
    int n;
    std::vector<double> matrix;

public:
    CholeskySolver(const std::vector<double>& input_matrix, int size) 
        : matrix(input_matrix), n(size) {}

    // Метод Холецкого с использованием BLAS
    bool choleskyDecomposition() {
        for (int i = 0; i < n; ++i) {
            // Вычисляем диагональный элемент L[i][i]
            double diag = matrix[i * n + i];
            
            // Вычитаем сумму квадратов элементов строки i до текущего столбца
            if (i > 0) {
                // Используем BLAS для вычисления скалярного произведения
                // ddot вычисляет x^T * y
                diag -= cblas_ddot(i, &matrix[i * n], 1, &matrix[i * n], 1);
            }
            
            if (diag <= 0.0) {
                std::cerr << "Матрица не положительно определена!" << std::endl;
                return false;
            }
            
            matrix[i * n + i] = std::sqrt(diag);
            
            // Обновляем элементы под диагональю в столбце i
            if (i < n - 1) {
                for (int j = i + 1; j < n; ++j) {
                    double element = matrix[j * n + i];
                    
                    // Вычитаем скалярное произведение строк i и j до текущего столбца
                    if (i > 0) {
                        element -= cblas_ddot(i, &matrix[j * n], 1, &matrix[i * n], 1);
                    }
                    
                    matrix[j * n + i] = element / matrix[i * n + i];
                }
            }
        }
        
        // Обнуляем верхний треугольник
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                matrix[i * n + j] = 0.0;
            }
        }
        
        return true;
    }

    // Умножение матриц с использованием BLAS (dgemm)
    std::vector<double> matrixMultiply(const std::vector<double>& B) const {
        std::vector<double> result(n * n, 0.0);
        
        // C = alpha * A * B + beta * C
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                   n, n, n, 1.0, matrix.data(), n, B.data(), n, 0.0, result.data(), n);
        
        return result;
    }

    // Проверка разложения: L * L^T должно равняться исходной матрице
    std::vector<double> verifyDecomposition(const std::vector<double>& original) const {
        // Вычисляем L * L^T
        std::vector<double> L_transposed = getTransposed();
        std::vector<double> reconstructed = matrixMultiply(L_transposed);
        
        // Вычисляем разницу с исходной матрицей
        std::vector<double> difference(n * n);
        for (int i = 0; i < n * n; ++i) {
            difference[i] = std::abs(reconstructed[i] - original[i]);
        }
        
        return difference;
    }

    // Получить транспонированную матрицу (L^T)
    std::vector<double> getTransposed() const {
        std::vector<double> transposed(n * n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                transposed[i * n + j] = matrix[j * n + i];
            }
        }
        return transposed;
    }

    // Вывод матрицы
    void printMatrix(const std::string& name = "") const {
        if (!name.empty()) {
            std::cout << name << ":\n";
        }
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                std::cout << matrix[i * n + j] << "\t";
            }
            std::cout << "\n";
        }
        std::cout << std::endl;
    }

    const std::vector<double>& getMatrix() const { return matrix; }
};

// Функция для создания случайной симметричной положительно определенной матрицы
std::vector<double> generateSPDMatrix(int n) {
    std::vector<double> A(n * n);
    std::vector<double> B(n * n);
    
    // Создаем случайную матрицу B
    for (int i = 0; i < n * n; ++i) {
        B[i] = static_cast<double>(rand()) / RAND_MAX;
    }
    
    // A = B * B^T + I (для обеспечения положительной определенности)
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
               n, n, n, 1.0, B.data(), n, B.data(), n, 0.0, A.data(), n);
    
    // Добавляем единичную матрицу для улучшения обусловленности
    for (int i = 0; i < n; ++i) {
        A[i * n + i] += 1.0;
    }
    
    return A;
}

int main() {
    const int n = 4;
    
    // Создаем симметричную положительно определенную матрицу
    std::vector<double> A = generateSPDMatrix(n);
    
    std::cout << "Однопараметрический метод Холецкого с использованием BLAS\n";
    std::cout << "=================================================\n\n";
    
    // Сохраняем исходную матрицу для проверки
    std::vector<double> original_A = A;
    
    // Создаем решатель и выполняем разложение Холецкого
    CholeskySolver solver(A, n);
    
    std::cout << "Исходная матрица A:\n";
    solver.printMatrix();
    
    if (solver.choleskyDecomposition()) {
        std::cout << "Нижнетреугольная матрица L после разложения Холецкого:\n";
        solver.printMatrix();
        
        // Проверка разложения
        std::vector<double> error = solver.verifyDecomposition(original_A);
        
        std::cout << "Максимальная ошибка восстановления: ";
        double max_error = 0.0;
        for (double err : error) {
            if (err > max_error) max_error = err;
        }
        std::cout << max_error << std::endl;
        
        if (max_error < 1e-10) {
            std::cout << "Разложение выполнено корректно!\n";
        } else {
            std::cout << "Внимание: возможная ошибка в разложении!\n";
        }
    }
    
    return 0;
}