#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdio>

using namespace std;

// Defines
#define PE_ARR_WIDTH = 16
#define ENABLE_DEBUG

float** fp_matrix;
int** i_matrix;
 
void check_inputs(char matrix_type, float** fp_matrix, int** i_matrix,
     int input_rows, int input_columns) {
    while(true){
        std::cout << "\nMatrix type:\n";
        printf("%c\n",matrix_type);
        std::cout << "\nMatrix contents:\n";
        if(matrix_type == 'f'){
        for (int i = 0; i < input_rows; ++i) {
            for (int j = 0; j < input_columns; ++j) {
                std::cout << fp_matrix[i][j] << " ";
            }
        std::cout << "\n";
        }
        } else if (matrix_type == 'i'){
            for (int i = 0; i < input_rows; ++i) {
                for (int j = 0; j < input_columns; ++j) {
                std::cout << i_matrix[i][j] << " ";
            }
        std::cout << "\n";            
        }
        }
        std::string input;
        std::cout << "Continue Printing?\n";
        std::getline(std::cin, input);
        if(input == "q"){
            break;
        }
    }
}

void cleanup_mem(){
// deallocate i_matrix
    if (i_matrix != nullptr) {
    for (int i = 0; i < input_rows; ++i) {
        delete[] i_matrix[i];
    }
    delete[] i_matrix;
    i_matrix = nullptr;
    }

// deallocate fp_matrix
if (fp_matrix != nullptr) {
    for (int i = 0; i < input_rows; ++i) {
        delete[] fp_matrix[i];
    }
    delete[] fp_matrix;
    fp_matrix = nullptr;
}
}

int main(){

    std::string filename;
    std::cout << "Enter the input file name: ";
    std::cin >> filename;

    std::ifstream infile(filename);
    if (!infile) {
        std::cerr << "Error: Cannot open file '" << filename << "'\n";
        return 1;
    }

    int input_rows, input_columns;
    char matrix_type;
    infile >> input_rows >> input_columns >> matrix_type;

    if(matrix_type == 'f'){
        fp_matrix = new float* [input_rows];
        for (int i = 0; i < input_rows; ++i) {
            fp_matrix[i] = new float[input_columns];
        }

    for (int i = 0; i < input_rows; ++i){
        for (int j = 0; j < input_columns; ++j){
            infile >> fp_matrix[i][j];
        }
    }
    infile.close();
    } else if(matrix_type == 'i') {
        i_matrix = new int* [input_rows];
        for (int i = 0; i < input_rows; ++i) {
            i_matrix[i] = new int[input_columns];
        }
            for (int i = 0; i < input_rows; ++i){
        for (int j = 0; j < input_columns; ++j){
            infile >> i_matrix[i][j];
        }
    }
    infile.close();
    } else {
        std::cerr << "Error: Invalid Matrix Type\n";
        return 1;
    }
    
    #ifdef ENABLE_DEBUG
    check_inputs(matrix_type,fp_matrix,i_matrix,input_rows,input_columns);
    #endif

    std::cout << "Exited Printing\n";

    cleanup_mem();

    std::cout << "Memory Freed. Exiting Program\n";
    return 0;
}