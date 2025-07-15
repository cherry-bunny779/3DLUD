#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdio>
#include <limits>
#include <iomanip>
#include "LU_decomp_blocked.h"

using namespace std;

// Defines
#define PE_ARR_WIDTH = 16
#define ENABLE_DEBUG
//#define VERBOSE

// DBG Defines
//#define DisableDBG_deep_copy_2d

// Global Vars
float **fp_matrix = nullptr;
int **i_matrix = nullptr;
int **perm_matrix = nullptr;
int **lmatrix_i = nullptr;
int **umatrix_i = nullptr;
float **lmatrix_fp = nullptr;
float **umatrix_fp = nullptr;
void* pivot_in = nullptr;
void* result_matrix = nullptr;

void print_matrix(char matrix_type, int **i_matrix, float **fp_matrix,
                  int input_rows, int input_columns)
{
#ifdef VERBOSE
    static int call = 0;
    std::cout << "print_matrix called: " << ++call << " times\n";
    std::cout << "\nMatrix type:\n";
    printf("%c\n", matrix_type);
#endif
    std::cout << "\nMatrix contents:\n";
    if (matrix_type == 'f')
    {
        for (int i = 0; i < input_rows; ++i)
        {
            for (int j = 0; j < input_columns; ++j)
            {
                std::cout << fp_matrix[i][j] << " ";
            }
            std::cout << "\n";
        }
    }
    else if (matrix_type == 'i')
    {
        for (int i = 0; i < input_rows; ++i)
        {
            for (int j = 0; j < input_columns; ++j)
            {
                std::cout << i_matrix[i][j] << " ";
            }
            std::cout << "\n";
        }
    }
}

void deep_copy_2d(int **src_i, float **src_f,
                  int **&dest_i, float **&dest_f,
                  int rows, int cols)
{
    if (src_i != nullptr) {
        dest_i = new int *[rows];
        for (int i = 0; i < rows; ++i) {
            dest_i[i] = new int[cols];
            for (int j = 0; j < cols; ++j) {
                dest_i[i][j] = src_i[i][j];
            }
        }
        dest_f = nullptr;
        std::cout << "Deep copied i_matrix.\n";
    }
    else if (src_f != nullptr) {
        dest_f = new float *[rows];
        for (int i = 0; i < rows; ++i) {
            dest_f[i] = new float[cols];
            for (int j = 0; j < cols; ++j) {
                dest_f[i][j] = src_f[i][j];
            }
        }
        dest_i = nullptr;
        std::cout << "Deep copied fp_matrix.\n";
    }
    else {
        std::cout << "Both source matrices are nullptr — nothing to copy.\n";
        dest_i = nullptr;
        dest_f = nullptr;
    }
}

void deallocate_2d(int **i_matrix_copy, float **fp_matrix_copy, int rows)
{
    if (i_matrix_copy != nullptr)
    {
        for (int i = 0; i < rows; ++i)
        {
            delete[] i_matrix_copy[i];
        }
        delete[] i_matrix_copy;
        std::cout << "Deallocated i_matrix_copy.\n";
    }

    if (fp_matrix_copy != nullptr)
    {
        for (int i = 0; i < rows; ++i)
        {
            delete[] fp_matrix_copy[i];
        }
        delete[] fp_matrix_copy;
        std::cout << "Deallocated fp_matrix_copy.\n";
    }

    if (i_matrix_copy == nullptr && fp_matrix_copy == nullptr)
    {
        std::cout << "Nothing to deallocate — both copies are nullptr.\n";
    }
}

// [To-Do][Optimize] Reuse deep_copy_2d(), pass i_matrix and fp_matrix instead
void partial_pivot(
    void* input_matrix, char matrix_type, int rows, int cols,
    void** result_matrix,              // output: pivoted matrix
    int** &permutation_matrix,        // output: permutation matrix
    int& row_operations               // output: number of row swaps
) {
    row_operations = 0;

    if (matrix_type == 'i') {
        int** in_mat = static_cast<int**>(input_matrix);
        int** res_mat = new int*[rows];
        int** perm_mat = new int*[rows];

        // Allocate and copy input matrix & identity permutation matrix
        for (int i = 0; i < rows; ++i) {
            res_mat[i] = new int[cols];
            perm_mat[i] = new int[rows];
            for (int j = 0; j < cols; ++j)
                res_mat[i][j] = in_mat[i][j];
            for (int j = 0; j < rows; ++j)
                perm_mat[i][j] = (i == j ? 1 : 0);
        }

        for (int k = 0; k < std::min(rows, cols); ++k) {
            int max_row = k;
            for (int i = k + 1; i < rows; ++i) {
                if (std::abs(res_mat[i][k]) > std::abs(res_mat[max_row][k])) {
                    max_row = i;
                }
            }
            if (max_row != k) {
                std::swap(res_mat[k], res_mat[max_row]);
                std::swap(perm_mat[k], perm_mat[max_row]);
                ++row_operations;
            }
        }

        *result_matrix = res_mat;
        permutation_matrix = perm_mat;

    } else if (matrix_type == 'f') {
        double** in_mat = static_cast<double**>(input_matrix);
        double** res_mat = new double*[rows];
        int** perm_mat = new int*[rows];

        // Allocate and copy input matrix & identity permutation matrix
        for (int i = 0; i < rows; ++i) {
            res_mat[i] = new double[cols];
            perm_mat[i] = new int[rows];
            for (int j = 0; j < cols; ++j)
                res_mat[i][j] = in_mat[i][j];
            for (int j = 0; j < rows; ++j)
                perm_mat[i][j] = (i == j ? 1.0 : 0.0);
        }

        for (int k = 0; k < std::min(rows, cols); ++k) {
            int max_row = k;
            for (int i = k + 1; i < rows; ++i) {
                if (std::abs(res_mat[i][k]) > std::abs(res_mat[max_row][k])) {
                    max_row = i;
                }
            }
            if (max_row != k) {
                std::swap(res_mat[k], res_mat[max_row]);
                std::swap(perm_mat[k], perm_mat[max_row]);
                ++row_operations;
            }
        }

        *result_matrix = res_mat;
        permutation_matrix = perm_mat;

    } else {
        std::cerr << "Unsupported matrix type '" << matrix_type << "'\n";
        *result_matrix = nullptr;
        permutation_matrix = nullptr;
    }
}

void lu_decomposition(
    int** src_i,
    float** src_f,
    char matrix_type,
    int n,
    int**& L_int,      // int version
    int**& U_int,
    float**& L_float,  // float version
    float**& U_float,
    int& num_mac,
    int& num_subt,
    int& num_div
) {
    num_mac = 0;
    num_subt = 0;
    num_div = 0;

    if (matrix_type == 'i') {
        int** A = src_i;
        L_int = new int*[n];
        U_int = new int*[n];

        for (int i = 0; i < n; ++i) {
            L_int[i] = new int[n];
            U_int[i] = new int[n];
            for (int j = 0; j < n; ++j) {
                L_int[i][j] = 0;
                U_int[i][j] = 0;
            }
        }

        for (int k = 0; k < n; ++k) {
            for (int j = k; j < n; ++j) {
                int sum = 0;
                for (int s = 0; s < k; ++s) {
                    sum += L_int[k][s] * U_int[s][j];
                    num_mac++;
                }
                U_int[k][j] = A[k][j] - sum;
                num_subt++;
            }

            L_int[k][k] = 1;

            for (int i = k + 1; i < n; ++i) {
                int sum = 0;
                for (int s = 0; s < k; ++s) {
                    sum += L_int[i][s] * U_int[s][k];
                    num_mac++;
                }
                if (U_int[k][k] == 0) {
                    std::cerr << "Zero pivot encountered at row " << k << "\n";
                    printf("Exiting the LUD Function\n");
                    return;
                }
                L_int[i][k] = (A[i][k] - sum) / U_int[k][k];
                num_div++;
                num_subt++;
            }
        }

    } else if (matrix_type == 'f') {
        float** A = src_f;
        L_float = new float*[n];
        U_float = new float*[n];

        for (int i = 0; i < n; ++i) {
            L_float[i] = new float[n];
            U_float[i] = new float[n];
            for (int j = 0; j < n; ++j) {
                L_float[i][j] = 0.0f;
                U_float[i][j] = 0.0f;
            }
        }

        for (int k = 0; k < n; ++k) {
            for (int j = k; j < n; ++j) {
                float sum = 0.0f;
                for (int s = 0; s < k; ++s) {
                    sum += L_float[k][s] * U_float[s][j];
                    num_mac++;
                }
                U_float[k][j] = A[k][j] - sum;
                num_subt++;
            }

            L_float[k][k] = 1.0f;

            for (int i = k + 1; i < n; ++i) {
                float sum = 0.0f;
                for (int s = 0; s < k; ++s) {
                    sum += L_float[i][s] * U_float[s][k];
                    num_mac++;
                }
                if (U_float[k][k] == 0.0f) {
                    std::cerr << "Zero pivot encountered at row " << k << "\n";
                    exit(1);
                }
                L_float[i][k] = (A[i][k] - sum) / U_float[k][k];
                num_div++;
                num_subt++;
            }
        }

    } else {
        std::cerr << "Invalid matrix type '" << matrix_type << "' (use 'i' or 'f')\n";
        exit(1);
    }
}

void cleanup_mem(int input_rows) {
    // deallocate i_matrix
    if (i_matrix != nullptr) {
        for (int i = 0; i < input_rows; ++i) {
            delete[] i_matrix[i];
        }
        delete[] i_matrix;
        i_matrix = nullptr;
        printf("Deallocated i_matrix.\n");
    }

    // deallocate fp_matrix
    if (fp_matrix != nullptr) {
        for (int i = 0; i < input_rows; ++i) {
            delete[] fp_matrix[i];
        }
        delete[] fp_matrix;
        fp_matrix = nullptr;
        printf("Deallocated fp_matrix.\n");
    }

    // deallocate perm_matrix
    if (perm_matrix != nullptr) {
        for (int i = 0; i < input_rows; ++i) {
            delete[] perm_matrix[i];
        }
        delete[] perm_matrix;
        perm_matrix = nullptr;
        printf("Deallocated perm_matrix.\n");
    }

    // deallocate lmatrix_i
    if (lmatrix_i != nullptr) {
        for (int i = 0; i < input_rows; ++i) {
            delete[] lmatrix_i[i];
        }
        delete[] lmatrix_i;
        lmatrix_i = nullptr;
        printf("Deallocated lmatrix_i.\n");
    }

    // deallocate umatrix_i
    if (umatrix_i != nullptr) {
        for (int i = 0; i < input_rows; ++i) {
            delete[] umatrix_i[i];
        }
        delete[] umatrix_i;
        umatrix_i = nullptr;
        printf("Deallocated umatrix_i.\n");
    }

    // deallocate lmatrix_fp
    if (lmatrix_fp != nullptr) {
        for (int i = 0; i < input_rows; ++i) {
            delete[] lmatrix_fp[i];
        }
        delete[] lmatrix_fp;
        lmatrix_fp = nullptr;
        printf("Deallocated lmatrix_fp.\n");
    }

    // deallocate umatrix_fp
    if (umatrix_fp != nullptr) {
        for (int i = 0; i < input_rows; ++i) {
            delete[] umatrix_fp[i];
        }
        delete[] umatrix_fp;
        umatrix_fp = nullptr;
        printf("Deallocated umatrix_fp.\n");
    }

    // deallocate pivot_in if allocated (void*)
    if (pivot_in != nullptr) {
        delete[] static_cast<char*>(pivot_in); // treat as block of memory
        pivot_in = nullptr;
        printf("Deallocated pivot_in.\n");
    }

    // deallocate result_matrix if allocated (void*)
    if (result_matrix != nullptr) {
        delete[] static_cast<char*>(result_matrix); // treat as block of memory
        result_matrix = nullptr;
        printf("Deallocated result_matrix.\n");
    }

    std::cout << "Memory Freed. Exiting Program\n";
}

void operation1()
{
    printf("Operation1\n");
}

void operation2()
{
    printf("Operation2\n");
}

void operation3()
{
    printf("Operation3\n");
}

int main()
{

    int input_rows, input_columns;
    char matrix_type;

    std::cout << "Welcome to 3DLUD-SIM!\n";

    while (true)
    {
        std::string filename;
        std::cout << "Enter the input file name: ";
        std::cin >> filename;

        if (filename == "quit")
        {
            std::cout << "Exiting 3DLUD-SIM\n";
            break;
        }
        std::ifstream infile(filename);
        if (!infile)
        {
            std::cerr << "Error: Cannot open file '" << filename << "'\n";
            continue;
        }

        infile >> input_rows >> input_columns >> matrix_type;

        if (matrix_type == 'f')
        {
            fp_matrix = new float *[input_rows];
            for (int i = 0; i < input_rows; ++i)
            {
                fp_matrix[i] = new float[input_columns];
            }

            for (int i = 0; i < input_rows; ++i)
            {
                for (int j = 0; j < input_columns; ++j)
                {
                    infile >> fp_matrix[i][j];
                }
            }
            infile.close();
        }
        else if (matrix_type == 'i')
        {
            i_matrix = new int *[input_rows];
            for (int i = 0; i < input_rows; ++i)
            {
                i_matrix[i] = new int[input_columns];
            }
            for (int i = 0; i < input_rows; ++i)
            {
                for (int j = 0; j < input_columns; ++j)
                {
                    infile >> i_matrix[i][j];
                }
            }
            infile.close();
        }
        else
        {
            std::cerr << "Error: Invalid Matrix Type\n";
            return 1;
        }

// Check Input Matrix Content Read From File
#ifdef ENABLE_DEBUG
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        while (true)
        {
            print_matrix(matrix_type, i_matrix, fp_matrix, input_rows, input_columns);
            std::string input;
            std::cout << "Check Input: Continue Printing? (Type 'c' to continue)\n";
            std::getline(std::cin, input);
            if (input == "c")
            {
                continue;
            }
            else
            {
                std::cout << "Exiting print_matrix()\n";
                break;
            }
        }
#endif

        while (true)
        {
            std::string choice;
            std::cout << "\nSelect an operation:\n";
            std::cout << " 1 - LUD (no pivoting)\n";
            std::cout << " 2 - Block LUD (no pivoting)\n";
            std::cout << " 3 - operation3\n";
            std::cout << " dbg\n";
            std::cout << "Type 'exit' to go back to file selection.\n";
            std::cout << "Your choice: ";

            std::getline(std::cin, choice);

            if (choice == "exit")
            {
                std::cout << "Returning to file selection.\n";
                break;
            }
            else if (choice == "1")
            {
                int num_row_ops,num_mac,num_subt,num_div;
                if(matrix_type == 'i'){ // this wrapper is needed for partial_pivot(). could use some optimization
                    pivot_in = static_cast<void*>(i_matrix);
                } else if(matrix_type == 'f'){
                    pivot_in = static_cast<void*>(fp_matrix);
                } else {
                    std::cerr << "Unsupported matrix type '" << matrix_type << "'\n";
                }
                partial_pivot(pivot_in,matrix_type,input_rows,input_columns,&result_matrix,perm_matrix,num_row_ops);
                #ifdef VERBOSE
                printf("Printing pivoted matrix\n");
                print_matrix(matrix_type,static_cast<int**>(result_matrix),static_cast<float**>(result_matrix),input_rows,input_columns);
                printf("\nPrinting Permutation Matrix\n");
                print_matrix('i',perm_matrix,nullptr,input_rows,input_columns);
                printf("\nv.s. Original Input\n");
                print_matrix(matrix_type,i_matrix,fp_matrix,input_rows,input_columns);
                #endif
                printf("\nPivoting: Number of Row Operations\n%i\n",num_row_ops);
                lu_decomposition(i_matrix,fp_matrix,matrix_type,input_rows,lmatrix_i,umatrix_i,lmatrix_fp,umatrix_fp,num_mac,num_subt,num_div);
                printf("\nNumber of Operations by Type:\nMAC: %i\nSubtraction: %i\nDivision: %i\n",num_mac,num_subt,num_div);
                if(matrix_type == 'i'){
                    printf("\nL Matrix\n");
                    print_matrix('i',lmatrix_i,lmatrix_fp,input_rows,input_columns);
                    printf("\nU Matrix\n");
                    print_matrix('i',umatrix_i,umatrix_fp,input_rows,input_columns);
                } else if(matrix_type == 'f'){
                    printf("\nL Matrix\n");
                    print_matrix('f',lmatrix_i,lmatrix_fp,input_rows,input_columns);
                    printf("\nU Matrix\n");
                    print_matrix('f',umatrix_i,umatrix_fp,input_rows,input_columns);
                } else {
                    std::cerr << "Unsupported matrix type '" << matrix_type << "'\n";
                }
                
            }
            else if (choice == "2")
            {
                operation2();
                float **fp_matrix_copy;
                int **i_matrix_copy; // 
                deep_copy_2d(nullptr,fp_matrix,i_matrix_copy,fp_matrix_copy,input_rows,input_columns);
                printf("\nComputing Blocked LUD\n");
                LU_Decomposition(fp_matrix_copy,input_rows,4);
                print_matrix('f',nullptr,fp_matrix_copy,input_rows,input_columns);
            }
            else if (choice == "3")
            {
                operation3();
            }
            else if (choice == "dbg")
            {
                #ifndef DisableDBG_deep_copy_2d
                printf("Debug deep_copy_2d\n");
                float **fp_matrix_copy;
                int **i_matrix_copy;
                deep_copy_2d(i_matrix, fp_matrix, i_matrix_copy, fp_matrix_copy, input_rows, input_columns);
                print_matrix(matrix_type, i_matrix_copy, fp_matrix_copy, input_rows, input_columns);
                deallocate_2d(i_matrix_copy, fp_matrix_copy, input_rows);
                #endif
            }
            else
            {
                std::cout << "Invalid choice. Please enter 1, 2, 3, or exit.\n";
            }
        }
    }
    cleanup_mem(input_rows);

    return 0;
}