#include <cstdio>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <lapacke.h>

void generateMatrix(double* A, int n) {
    double* B = new double[n * n];
    for (int i = 0; i < n * n; i++) {
        B[i] = (double)(rand() % 100) / 10.0;
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int k = 0; k < n; k++) {
                sum += B[k * n + i] * B[k * n + j];
            }
            A[i * n + j] = sum;
        }
    }
    delete[] B;
}

void printMatrix(const char* name, double* A, int n) {
    printf("\n%s:\n", name);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%10.4f ", A[i * n + j]);
        }
        printf("\n");
    }
}

void printVector(const char* name, double* v, int n) {
    printf("\n%s:\n", name);
    for (int i = 0; i < n; i++) {
        printf("%10.6f\n", v[i]);
    }
}

void matrixVectorMult(double* A, double* v, double* result, int n) {
    for (int i = 0; i < n; i++) {
        result[i] = 0.0;
        for (int j = 0; j < n; j++) {
            result[i] += A[i * n + j] * v[j];
        }
    }
}

double vectorNorm(double* v, int n) {
    double norm = 0.0;
    for (int i = 0; i < n; i++) {
        norm += v[i] * v[i];
    }
    return sqrt(norm);
}

void choleskyDecomposition(double* A, int n) {
    printf("\nРазложение Холецкого\n");

    double* A_copy = new double[n * n];
    memcpy(A_copy, A, n * n * sizeof(double));

    lapack_int info = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'L', n, A_copy, n);

    if (info == 0) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (j <= i) {
                    printf("%10.4f ", A_copy[i * n + j]);
                } else {
                    printf("%10.4f ", 0.0);
                }
            }
            printf("\n");
        }
    }
    delete[] A_copy;
}

void svdDecomposition(double* A, int n) {
    printf("\nSVD разложение\n");

    double* A_copy = new double[n * n];
    memcpy(A_copy, A, n * n * sizeof(double));

    double* U = new double[n * n];
    double* S = new double[n];
    double* VT = new double[n * n];

    lapack_int info;
    double* superb = new double[n - 1];
    info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'A', 'A', n, n, A_copy, n, S, U, n, VT, n, superb);

    if (info == 0) {
        printf("S:\n");
        for (int i = 0; i < n; i++) {
            printf("%10.6f\n", S[i]);
        }

        printf("\nU:\n");
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                printf("%10.4f ", U[i * n + j]);
            }
            printf("\n");
        }

        printf("\nV^T:\n");
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                printf("%10.4f ", VT[i * n + j]);
            }
            printf("\n");
        }
    }

    delete[] A_copy;
    delete[] U;
    delete[] S;
    delete[] VT;
    delete[] superb;
}

void qrDecomposition(double* A, int n) {
    printf("\nQR разложение\n");

    double* A_copy = new double[n * n];
    memcpy(A_copy, A, n * n * sizeof(double));

    double* tau = new double[n];
    lapack_int info = LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, n, n, A_copy, n, tau);

    if (info == 0) {
        printf("R:\n");
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (j >= i) {
                    printf("%10.4f ", A_copy[i * n + j]);
                } else {
                    printf("%10.4f ", 0.0);
                }
            }
            printf("\n");
        }

        double* Q = new double[n * n];
        memcpy(Q, A_copy, n * n * sizeof(double));

        info = LAPACKE_dorgqr(LAPACK_ROW_MAJOR, n, n, n, Q, n, tau);

        if (info == 0) {
            printf("\nQ:\n");
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    printf("%10.4f ", Q[i * n + j]);
                }
                printf("\n");
            }
        }
        delete[] Q;
    }

    delete[] A_copy;
    delete[] tau;
}

void eigenvaluesDecomposition(double* A, int n) {
    printf("\nСобственные значения и векторы\n");

    double* A_copy = new double[n * n];
    memcpy(A_copy, A, n * n * sizeof(double));

    double* W = new double[n];
    lapack_int info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', n, A_copy, n, W);

    if (info == 0) {
        printf("λ:\n");
        for (int i = 0; i < n; i++) {
            printf("%10.6f\n", W[i]);
        }

        printf("\nV:\n");
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                printf("%10.6f ", A_copy[j * n + i]);
            }
            printf("\n");
        }

        double* v = new double[n];
        double* Av = new double[n];

        for (int i = 0; i < n; i++) {
            v[i] = A_copy[i * n + 0];
        }

        double* A_original = new double[n * n];
        memcpy(A_original, A, n * n * sizeof(double));

        matrixVectorMult(A_original, v, Av, n);

        printf("\nA*v:\n");
        for (int i = 0; i < n; i++) {
            printf("%10.6f ", Av[i]);
        }
        printf("\n");

        printf("\nλ*v:\n");
        for (int i = 0; i < n; i++) {
            printf("%10.6f ", W[0] * v[i]);
        }
        printf("\n");

        delete[] v;
        delete[] Av;
        delete[] A_original;
    }
    delete[] A_copy;
    delete[] W;
}

void invertMatrix(double* A, int n) {
    printf("\nОбращение матрицы\n");

    double* A_copy = new double[n * n];
    memcpy(A_copy, A, n * n * sizeof(double));

    lapack_int* ipiv = new lapack_int[n];
    lapack_int info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, A_copy, n, ipiv);

    if (info == 0) {
        info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, A_copy, n, ipiv);
        if (info == 0) {
            printMatrix("A^{-1}", A_copy, n);
        } else {
            printf("Ошибка dgetri: %d\n", info);
        }
    }

    delete[] A_copy;
    delete[] ipiv;
}

int main() {
    int n = 4;

    double* A = new double[n * n];
    generateMatrix(A, n);

    printMatrix("A", A, n);

    double* b = new double[n];
    for (int i = 0; i < n; i++) {
        b[i] = (double)(rand() % 100) / 10.0;
    }
    printVector("b", b, n);

    choleskyDecomposition(A, n);
    svdDecomposition(A, n);
    qrDecomposition(A, n);
    eigenvaluesDecomposition(A, n);
    invertMatrix(A, n);

    delete[] A;
    delete[] b;

    return 0;
}
