#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdio>
#include <limits>

using namespace std;

// Defines
#define PE_ARR_WIDTH = 16
#define ENABLE_DEBUG

float** fp_matrix;
int** i_matrix;
 
void check_inputs(char matrix_type, float** fp_matrix, int** i_matrix,
     int input_rows, int input_columns) {

        static int call = 0;
        std::cout << "check_inputs called: " << ++call << " times\n";
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
            printf("Printing i_matrix\n");
            for (int i = 0; i < input_rows; ++i) {
                for (int j = 0; j < input_columns; ++j) {
                std::cout << i_matrix[i][j] << " ";
            }
        std::cout << "\n";            
        }
        }
    
}

void cleanup_mem(int input_rows){
// deallocate i_matrix
    if (i_matrix != nullptr) {
    for (int i = 0; i < input_rows; ++i) {
        delete[] i_matrix[i];
    }
    delete[] i_matrix;
    i_matrix = nullptr;
    #ifdef ENABLE_DEBUG
        printf("Deallocated i_matrix.\n");
    #endif
    }

// deallocate fp_matrix
if (fp_matrix != nullptr) {
    for (int i = 0; i < input_rows; ++i) {
        delete[] fp_matrix[i];
    }
    delete[] fp_matrix;
    fp_matrix = nullptr;
    #ifdef ENABLE_DEBUG
        printf("Deallocated fp_matrix.\n");
    #endif
}
    std::cout << "Memory Freed. Exiting Program\n";
}

void operation1(){
    printf("Operation1\n");
}

void operation2(){
    printf("Operation2\n");
}

void operation3(){
    printf("Operation3\n");
}

int main(){

    int input_rows, input_columns;
    char matrix_type;

    std::cout << "Welcome to 3DLUD-SIM!\n";

    while(true){
    std::string filename;
    std::cout << "Enter the input file name: ";
    std::cin >> filename;

    if(filename == "quit"){
        std::cout << "Exiting 3DLUD-SIM\n";
        break;
    }
    std::ifstream infile(filename);
    if (!infile) {
        std::cerr << "Error: Cannot open file '" << filename << "'\n";
        continue;
    }

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
    
    // Check Input Matrix Content Read From File
    #ifdef ENABLE_DEBUG
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    while(true){
        check_inputs(matrix_type,fp_matrix,i_matrix,input_rows,input_columns);
        std::string input;
        std::cout << "Check Input: Continue Printing? (Type 'c' to continue)\n";
        std::getline(std::cin, input);
        if(input == "c"){
            continue;
        } else {
            std::cout << "Exiting check_inputs()\n";
            break;
        }
    }
    #endif
            while (true) {
            std::string choice;
            std::cout << "\nSelect an operation:\n";
            std::cout << " 1 - operation1\n";
            std::cout << " 2 - operation2\n";
            std::cout << " 3 - operation3\n";
            std::cout << "Type 'exit' to go back to file selection.\n";
            std::cout << "Your choice: ";

            std::getline(std::cin, choice);

            if (choice == "exit") {
                std::cout << "Returning to file selection.\n";
                break;
            } else if (choice == "1") {
                operation1();
            } else if (choice == "2") {
                operation2();
            } else if (choice == "3") {
                operation3();
            } else {
                std::cout << "Invalid choice. Please enter 1, 2, 3, or exit.\n";
            }
        }
}
    cleanup_mem(input_rows);

    return 0;
}