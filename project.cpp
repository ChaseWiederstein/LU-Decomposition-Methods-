// FCM_HW2.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <chrono>

std::vector<std::vector<double>> NoPivot(std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& L) {
    int n = A.size();
    L.resize(n, std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> U(n, std::vector<double>(n, 0.0));

    for (int i = 0; i < n; i++) {
        L[i][i] = 1;
        for (int k = i; k < n; k++) {
            U[i][k] = A[i][k];
            for (int j = 0; j < i; j++) {
                U[i][k] -= L[i][j] * U[j][k];
            }
        }
        for (int k = i + 1; k < n; k++) {
            L[k][i] = A[k][i];
            for (int j = 0; j < i; j++) {
                L[k][i] -= L[k][j] * U[j][i];
            }
            L[k][i] /= U[i][i];
        }
    }

    return U;
}

std::vector<std::vector<double>> PartialPivot(std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& L, std::vector<std::vector<double>>& P, std::vector<double>& pivoting_vector) {
    int n = A.size();
    L.resize(n, std::vector<double>(n, 0.0));
    P = std::vector<std::vector<double>>(n, std::vector<double>(n, 0));
    pivoting_vector.resize(n);

    // Initialize P as the identity matrix
    for (int i = 0; i < n; i++) {
        P[i][i] = 1;
        pivoting_vector[i] = i;
    }

    std::vector<std::vector<double>> U = A;  // Initialize U as a copy of A (identity-like)

    for (int k = 0; k < n - 1; k++) {
        int max_row = k;

        // Partial Pivoting: Find the row with the largest absolute value in the current column
        for (int i = k; i < n; i++) {
            if (std::abs(U[i][k]) > std::abs(U[max_row][k])) {
                max_row = i;
            }
        }

        // Swap rows in U, P, and the pivoting vector
        if (max_row != k) {
            std::swap(U[k], U[max_row]);
            std::swap(P[k], P[max_row]);
            std::swap(pivoting_vector[k], pivoting_vector[max_row]);
        }

        // Perform LU decomposition
        for (int i = k + 1; i < n; i++) {
            L[i][k] = U[i][k] / U[k][k];
            for (int j = k; j < n; j++) {
                U[i][j] -= L[i][k] * U[k][j];
            }
        }
    }

    // Populate the diagonal of L with 1s
    for (int i = 0; i < n; i++) {
        L[i][i] = 1.0;
    }

    return U; // Return the upper triangular matrix U
}

void CompletePivot(std::vector<std::vector<int>>& A, std::vector<std::vector<int>>& L, std::vector<std::vector<int>>& U, std::vector<std::vector<int>>& P, std::vector<std::vector<int>>& Q) {
    int n = A.size();

    // Initialize the permutation matrices P and Q as identity matrices
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            P[i][j] = Q[i][j] = (i == j) ? 1 : 0;
        }
    }

    for (int i = 0; i < n; i++) {
        int maxVal = 0;
        int pivotRow = -1;
        int pivotCol = -1;

        // Find the pivot element with complete pivoting
        for (int k = i; k < n; k++) {
            for (int l = i; l < n; l++) {
                if (abs(A[k][l]) > maxVal) {
                    maxVal = abs(A[k][l]);
                    pivotRow = k;
                    pivotCol = l;
                }
            }
        }

        // Swap rows in A and P
        std::swap(A[i], A[pivotRow]);
        std::swap(P[i], P[pivotRow]);

        // Swap columns in A and Q
        for (int row = 0; row < n; row++) {
            std::swap(A[row][i], A[row][pivotCol]);
            std::swap(Q[row][i], Q[row][pivotCol]);
        }

        L[i][i] = 1;
        for (int j = i; j < n; j++) {
            U[i][j] = A[i][j];
            for (int k = 0; k < i; k++) {
                U[i][j] -= L[i][k] * U[k][j];
            }
        }

        for (int j = i + 1; j < n; j++) {
            L[j][i] = A[j][i];
            for (int k = 0; k < i; k++) {
                L[j][i] -= L[j][k] * U[k][i];
            }
            L[j][i] /= U[i][i];
        }
    }
}

std::vector<double> forwardSubstitution(const std::vector<std::vector<double>>& L, const std::vector<double>& b) {
    int n = L.size();
    std::vector<double> y(n, 0);

    for (int i = 0; i < n; i++) {
        if (L[i][i] == 0) {
            std::cout << "ERROR: Det(A) is nonzero" << std::endl;
        }
        else {
            y[i] = b[i];
            for (int j = 0; j < i; j++) {
                y[i] -= L[i][j] * y[j];
            }
            y[i] /= L[i][i];
        }
    }

    return y;
}

std::vector<double> backwardSubstitution(const std::vector<std::vector<double>>& U, const std::vector<double>& y) {
    int n = U.size();
    std::vector<double> x_tilde(n, 0);

    for (int i = n - 1; i >= 0; i--) {
        if (U[i][i] == 0) {
            std::cout << "ERROR: Det(A) is nonzero" << std::endl;
        }
        else {
            x_tilde[i] = y[i];
            for (int j = n - 1; j > i; j--) {
                x_tilde[i] -= U[i][j] * x_tilde[j];
            }
            x_tilde[i] /= U[i][i];
        }

    }

    return x_tilde;
}

std::vector<std::vector<double>> matrixMatrixMultiplication(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B) {
    int numRowsA = A.size();
    int numColsA = A[0].size();
    int numRowsB = B.size();
    int numColsB = B[0].size();

    if (numColsA != numRowsB) {
        std::cout << "Matrix dimensions are incompatible for multiplication." << std::endl;
        return std::vector<std::vector<double>>();
    }

    std::vector<std::vector<double>> result(numRowsA, std::vector<double>(numColsB, 0));

    for (int i = 0; i < numRowsA; i++) {
        for (int j = 0; j < numColsB; j++) {
            for (int k = 0; k < numColsA; k++) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return result;
}

std::vector<std::vector<double>> matrixSubtraction(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B) {
    int numRowsA = A.size();
    int numColsA = A[0].size();
    int numRowsB = B.size();
    int numColsB = B[0].size();

    if (numRowsA != numRowsB || numColsA != numColsB) {
        std::cout << "Matrix dimensions are incompatible for subtraction." << std::endl;
        return std::vector<std::vector<double>>();
    }

    std::vector<std::vector<double>> result(numRowsA, std::vector<double>(numColsA, 0));

    for (int i = 0; i < numRowsA; i++) {
        for (int j = 0; j < numColsA; j++) {
            result[i][j] = A[i][j] - B[i][j];
        }
    }

    return result;
}

std::vector<double> subtractVectors(const std::vector<double>& a, const std::vector<double>& b) {
    int size = a.size();
    std::vector<double> result(size);

    if (size != b.size()) {
        // Handle error or return an appropriate value
        return result; // A vector of zeros or empty, depending on your error-handling approach
    }

    for (int i = 0; i < size; i++) {
        result[i] = a[i] - b[i];
    }

    return result;
}

double matrixOneNorm(const std::vector<std::vector<double>>& A) {
    int numRows = A.size();
    int numCols = A[0].size();  // Assuming all columns have the same size

    double maxColumnSum = 0.0;

    for (int j = 0; j < numCols; j++) {
        double columnSum = 0.0;
        for (int i = 0; i < numRows; i++) {
            columnSum += std::abs(A[i][j]);
        }
        maxColumnSum = std::max(maxColumnSum, columnSum);
    }

    return maxColumnSum;
}

double generateRandomDouble(double min, double max) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<double> distribution(min, max);
    return distribution(gen);
}

std::vector<double> matrixVectorMultiply(const std::vector<std::vector<double>>& A, const std::vector<double>& x) {
    int numRows = A.size();
    int numCols = A[0].size();  // Assuming all columns have the same size

    std::vector<double> Ax(numRows, 0);

    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            Ax[i] += A[i][j] * x[j];
        }
    }

    return Ax;
}

std::vector<double> generateRandomDoubleVector(int length, double min, double max) {
    std::vector<double> vec(length);
    for (int i = 0; i < length; i++) {
        vec[i] = generateRandomDouble(min, max);
    }
    return vec;
}

double frobeniusNorm(const std::vector<std::vector<double>>& A) {
    double norm = 0.0;
    for (const auto& row : A) {
        for (int element : row) {
            norm += std::pow(element, 2);
        }
    }
    return std::sqrt(norm);
}

double vectorOneNorm(const std::vector<double>& x) {
    double norm = 0.0;
    for (int i = 0; i < x.size(); i++) {
        norm += std::abs(x[i]);
    }
    return norm;
}

double vectorTwoNorm(const std::vector<double>& x) {
    double norm = 0.0;
    for (int i = 0; i < x.size(); i++) {
        norm += std::pow(x[i], 2);
    }
    return std::sqrt(norm);
}

//----------------- SPD ---------------- 
std::vector<std::vector<double>> generateLowerTriangularMatrix(int n, double min, double max) {
    std::vector<std::vector<double>> matrix(n, std::vector<double>(n, 0.0));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            matrix[i][j] = generateRandomDouble(min, max);
        }
    }

    return matrix;
}

std::vector<std::vector<double>> transposeMatrix(const std::vector<std::vector<double>>& matrix) {
    int numRows = matrix.size();
    int numCols = matrix[0].size();
    std::vector<std::vector<double>> transpose(numCols, std::vector<double>(numRows, 0.0));

    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            transpose[j][i] = matrix[i][j];
        }
    }

    return transpose;
}

//-----------Test A for diagonaly dominant ---------------
bool isDiagonallyDominant(const std::vector<std::vector<double>>& A) {
    int n = A.size();

    for (int i = 0; i < n; i++) {
        double diagonal = std::abs(A[i][i]);
        double rowSum = 0.0;

        for (int j = 0; j < n; j++) {
            if (j != i) {
                rowSum += std::abs(A[i][j]);
            }
        }

        if (diagonal <= rowSum) {
            return false;
        }
    }

    return true;
}

//------------ Make a Diagonally Dominant Matrix A -------------------
std::vector<std::vector<double>> generateDiagonallyDominantMatrix(int n) {
    std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                // Set diagonal elements to random positive values
                A[i][j] = generateRandomDouble(1.0, 10.0); // Adjust the range as needed
            }
            else {
                // Set off-diagonal elements to random values
                A[i][j] = generateRandomDouble(-1.0, 1.0); // Adjust the range as needed
            }
        }
    }

    // Make the diagonal elements dominant
    for (int i = 0; i < n; i++) {
        double rowSum = 0.0;
        for (int j = 0; j < n; j++) {
            if (i != j) {
                rowSum += std::abs(A[i][j]);
            }
        }
        A[i][i] = 2.0 * rowSum; // You can adjust the factor as needed
    }

    return A;
}

std::vector<std::vector<double>> generateWellConditionedMatrix(int n, double conditionNumber) {
    std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> distribution(-1.0, 1.0); // Adjust the range as needed

    // Generate a random matrix with values between -1 and 1
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = distribution(gen);
        }
    }

    // Set the diagonal elements to ensure positive definiteness
    for (int i = 0; i < n; i++) {
        A[i][i] = conditionNumber; // Adjust the condition number as needed
    }

    return A;
}

std::vector<std::vector<double>> generateWellConditionedL(int n) {
    std::vector<std::vector<double>> L(n, std::vector<double>(n, 0.0));

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> distribution(1.0, 10.0); // Adjust the range as needed

    // Generate diagonal elements of L
    for (int i = 0; i < n; i++) {
        L[i][i] = distribution(gen); // Ensure diagonal elements are not too small or too large
    }

    // Generate off-diagonal elements of L
    for (int i = 1; i < n; i++) {
        for (int j = 0; j < i; j++) {
            L[i][j] = distribution(gen); // Adjust the range as needed
        }
    }

    return L;
}


int main() {

    int n = 100;

    //-----------Random A-------------------

    /*std::vector<std::vector<double>> A(n, std::vector<double>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = generateRandomDouble(1.0, 20.0);
        }
    }*/


    //-------------- SPD -------------------


    //--------------------------------------------
    double min = 1;
    double max = 20;
    std::vector<std::vector<double>> L_tri = generateLowerTriangularMatrix(n, min, max);
    std::vector<std::vector<double>> L_tri_trans = transposeMatrix(L_tri);
    std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));
    //Below does not make it diagonally dominant for partial
   /* for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k <= i; k++) {
                A[i][j] += L_tri[i][k] * L_tri_trans[k][j];
            }
        }
    }*/
   // Below makes A diag dom for noPivot()
    for (int i = 0; i < n; i++) {
        double rowSum = 0.0;
        for (int j = 0; j < n; j++) {
            if (i != j) {
                // Adjust these values as needed
                A[i][j] = generateRandomDouble(-1.0, 1.0);
                rowSum += std::abs(A[i][j]);
            }
        }
        // Set the diagonal element to a value greater than the sum of others
        A[i][i] = rowSum + 1.0;
    }
   // -------------------------------------- -

        //--------- generates well conditioned A for partial ------------------------
        //double targetConditionNumber = 10.0; 
        //std::vector<std::vector<double>> A = generateWellConditionedMatrix(n, targetConditionNumber);
        // ---------------------------------------------------------------------------------
        /*std::vector<std::vector<double>> L_tri = generateWellConditionedL(n);
        std::vector<std::vector<double>> L_tri_trans = transposeMatrix(L_tri);
        std::vector<std::vector<double>> A = matrixMatrixMultiplication(L_tri, L_tri_trans);
         */

         // -------------- Make A Random and diagonally dominant --------------
         //std::vector<std::vector<double>> A = generateDiagonallyDominantMatrix(n);

         //------------ Generate wellconditione L for partial pivot() --------------
        // std::vector<std::vector<double>> A = { {2, 1, 0},{-4, 0, 4},{2, 5, -10} };
        //std::vector<double> b = { 3, 0, 17 };
        //std::vector<double> x;

    std::vector<double> x = generateRandomDoubleVector(n, 1.0, 20.0);
    std::vector<std::vector<double>> L;
    std::vector<std::vector<double>> U;
    std::vector<std::vector<double>> P;
    std::vector<double> pivoting_vector;
    std::vector<double> b = matrixVectorMultiply(A, x);


    if (isDiagonallyDominant(A)) {
        std::cout << "Matrix A is diagonally dominant so were using NoPivot()" << std::endl;

        std::vector<std::vector<double>> U = NoPivot(A, L);
        auto start = std::chrono::high_resolution_clock::now();
        std::vector<double> y = forwardSubstitution(L, b);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        std::cout << "Time taken for forward substitution: " << duration.count() << " microseconds" << std::endl;
        start = std::chrono::high_resolution_clock::now();
        std::vector<double> x_tilde = backwardSubstitution(U, y);
        end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

        std::cout << "Time taken for backward substitution: " << duration.count() << " microseconds" << std::endl;
        //----------- part 1 -------------
        std::vector<std::vector<double>> LU = matrixMatrixMultiplication(L, U);
        std::vector<std::vector<double>> A_minus_LU = matrixSubtraction(A, LU);
        double oneNormA_minus_LU = matrixOneNorm(A_minus_LU);
        double oneNormA = matrixOneNorm(A);
        std::cout << "one norm ||A-LU||/||A|| =" << oneNormA_minus_LU / oneNormA << std::endl;

        double frobNormA_minus_LU = frobeniusNorm(A_minus_LU);
        double frobNormA = frobeniusNorm(A);
        std::cout << "Frob norm ||A-LU||/||A|| =" << frobNormA_minus_LU / frobNormA << std::endl;

        //------------- Part 2 --------------------------
        std::vector<double> x_minus_xtilde = subtractVectors(x, x_tilde);
        double x_minus_xtilde_oneNorm = vectorOneNorm(x_minus_xtilde);
        double xOneNorm = vectorOneNorm(x);
        std::cout << "One norm x: " << xOneNorm << std::endl;
        double relativeErrorOneNorm = (x_minus_xtilde_oneNorm) / xOneNorm;
        std::cout << "One Norm ||x-x~||/||x|| = " << relativeErrorOneNorm << std::endl;

        double x_minus_xtilde_TwoNorm = vectorTwoNorm(x_minus_xtilde);
        double xTwoNorm = vectorTwoNorm(x_tilde);
        double relativeErrorTwoNorm = x_minus_xtilde_TwoNorm / xTwoNorm;
        std::cout << "Two Norm ||x-x~||/||x|| = " << relativeErrorTwoNorm << std::endl;

        std::vector<double> A_x = matrixVectorMultiply(A, x_tilde);
        std::vector<double> b_minus_ax = subtractVectors(b, A_x);
        double b_minus_ax_oneNorm = vectorOneNorm(b_minus_ax);
        double bOneNorm = vectorOneNorm(b);
        double partThreeOneNorm = b_minus_ax_oneNorm / bOneNorm;
        std::cout << "One Norm ||b-Ax||/||b|| = " << partThreeOneNorm << std::endl;

        double b_minus_ax_TwoNorm = vectorTwoNorm(b_minus_ax);
        double bTwoNorm = vectorTwoNorm(b);
        double partThreeTwoNorm = b_minus_ax_TwoNorm / bTwoNorm;
        std::cout << "Two Norm ||b-Ax||/||b|| =" << partThreeTwoNorm << std::endl;


    }
    else {
        std::cout << "Matrix A is not diagonally dominant so were using PartialPivot()" << std::endl;
        std::vector<std::vector<double>> U = PartialPivot(A, L, P, pivoting_vector);
        // std::vector<double> b = { 8, 2, 3 };
        std::vector<double> y = forwardSubstitution(L, b);
        std::vector<double> x_tilde = backwardSubstitution(U, y);

        //---------- Part 1 ---------------------------
        std::vector<std::vector<double>> PA = matrixMatrixMultiplication(P, A);
        std::vector<std::vector<double>> LU = matrixMatrixMultiplication(L, U);
        std::vector<std::vector<double>> PA_minus_LU = matrixSubtraction(PA, LU);

        double oneNormPA_minus_LU = matrixOneNorm(PA_minus_LU);
        double oneNormA = matrixOneNorm(A);
        std::cout << "one norm ||PA-LU||/||A|| =" << oneNormPA_minus_LU / oneNormA << std::endl;
        double frobNormPA_minus_LU = frobeniusNorm(PA_minus_LU);
        double frobNormA = frobeniusNorm(A);
        std::cout << "Frob norm ||PALU||/||A||:" << frobNormPA_minus_LU / frobNormA << std::endl;

        //------------- Part 2 --------------------------
        std::vector<double> x_minus_xtilde = subtractVectors(x, x_tilde);
        double x_minus_xtilde_oneNorm = vectorOneNorm(x_minus_xtilde);
        double xOneNorm = vectorOneNorm(x);
        double relativeErrorOneNorm = x_minus_xtilde_oneNorm / xOneNorm;
        std::cout << "One-Norm ||x-x~||/||x|| = " << relativeErrorOneNorm << std::endl;

        double x_minus_xtilde_twoNorm = vectorTwoNorm(x_minus_xtilde);
        double xTwoNorm = vectorTwoNorm(x);
        double relativeErrorTwoNorm = x_minus_xtilde_twoNorm / xTwoNorm;
        std::cout << "Two Norm ||x-x~||/||x|| = " << relativeErrorTwoNorm << std::endl;

        std::vector<double> A_x = matrixVectorMultiply(A, x_tilde);
        std::vector<double> b_minus_Ax = subtractVectors(b, A_x);
        double b_minus_ax_oneNorm = vectorOneNorm(b_minus_Ax);
        double bOneNorm = vectorOneNorm(b);
        double partThreeOneNorm = b_minus_ax_oneNorm / bOneNorm;
        std::cout << "One norm ||b-Ax~||/||b|| =" << partThreeOneNorm << std::endl;

        double b_minus_ax_twoNorm = vectorTwoNorm(b_minus_Ax);
        double bTwoNorm = vectorTwoNorm(b);
        double partThreeTwoNorm = b_minus_ax_twoNorm / bTwoNorm;
        std::cout << "Two norm ||b-Ax~||/||b|| =" << partThreeTwoNorm << std::endl;

    }

    return 0;

}
