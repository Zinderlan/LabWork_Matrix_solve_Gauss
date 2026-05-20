#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <exception>
using std::vector;

//шаблонная функция матрицы
template <typename T>
class Matrix : public vector<vector<T>>
{
private:
    int count_strings; //количество строк в матрице
    int count_variable; //количество переменных(столбцов)

public:
    //конструктор инициализации
    Matrix(const int count_s, const int count_v) : vector<vector<T>>(count_s, std::vector<T>(count_v, T{ 0 })),
        count_strings{ count_s }, count_variable{ count_v }   { }

    T& operator()(int i, int j) {
        return (*this)[i][j];
    }
    const T& operator()(int i, int j) const {
        return (*this)[i][j];
    }

    //функция для введения строки
    void input_string(const int row_index)
    {
        vector<T> result(count_variable);
        std::string matrix_string;
        std::getline(std::cin >> std::ws, matrix_string);
        double temp;
        std::istringstream iss(matrix_string);
        for (int i = 0; i < count_variable && iss >> temp; i++)
        {
            T res = static_cast<T>(temp);
            result[i] = res;
        }
        (*this)[row_index] = std::move(result);
    }

    //функция печать матрицы
    void print_matrix()
    {
        std::cout << "Matrix: \n";
        for (int i = 0; i < count_strings; i++) {
            for (int j = 0; j < count_variable; j++) {
                std::cout << (*this)(i, j) << "\t";
            }
            std::cout << "\n";
        }
    }

    //функция поиска максимальноого элемента по col-му столбцу, начиная с start_row
    int Find_index_max(int start_row, int col)
    {
        T main_el{ (*this)(start_row, col) };
        int index_main{ start_row };

        for (int j = start_row; j < count_strings; j++)
        {
            //проверка максимального элемента
            if (std::abs(main_el) < std::abs((*this)(j, col)))
            {
                main_el = (*this)(j, col);
                index_main = j;
            }
        }
        return index_main;
    }

    //решение системы методом гаусса
    void solve_Gauss(vector<double>& free_coef, const vector<std::string>& Variable)
    {
        vector<int> lead_col_of_row(count_strings, -1);
        vector<int> free_variables(count_variable, 0);
        int free_count = 0;
        T epsilon = 1e-10;
        int row = 0;
        for (int col = 0; col < count_variable; col++)
        {
            int index_main = Find_index_max(row, col);
            T main_element{ (*this)(index_main, col) };

            if (std::abs(main_element) < epsilon) {
                continue;            //столбец пропускаем, row не увеличиваем
            }

            //swap между векторами, если главный не на row месте
            else if (index_main != row)
            {
                std::swap((*this)[row], (*this)[index_main]);
                std::swap(free_coef[row], free_coef[index_main]);
            }

            lead_col_of_row[row] = col;

            //деление ведущей строки(индекса row) и свободного члена
            free_coef[row] = free_coef[row] / main_element;
            for (int j = 0; j < count_variable; j++)
            {
                (*this)(row, j) = (*this)(row, j) / main_element;
            }

            //деление всех остальных строк и свобожгых членов, исключая row
            for (int k = 0; k < count_strings; k++)
            {
                if (k != row)
                {
                    T element = (*this)(k, col);
                    for (int j = 0; j < count_variable; j++)
                    {
                        (*this)(k, j) -= element * (*this)(row, j);
                    }
                    free_coef[k] -= element * free_coef[row];
                }
            }

            row++; //переход к следующей строке

        }

        //проверка на несоместимость системы
        for (int i = 0; i < count_strings; i++)
        {
            bool null_string = true;
            for (int j = 0; j < count_variable; j++)
            {
                if (std::abs((*this)(i, j)) > epsilon)
                {
                    null_string = false;
                    break;
                }
            }
            if (null_string && std::abs(free_coef[i]) > epsilon)
            {
                throw std::runtime_error("Error: inconsistent system of equations");
            }
        }

        //определяем свободные переменные: столбцы, не попавшие в lead_col_of_row
        for (int col = 0; col < count_variable; col++)
        {
            bool is_leading = false;
            for (int r = 0; r < row; r++)
            {
                if (lead_col_of_row[r] == col)
                {
                    is_leading = true;
                    break;
                }
            }
            if (not is_leading)
            {
                free_variables[col] = 1;   // помечаем как свободную
                free_count++;
            }
        }
        bool has_free = (free_count > 0);

        //вывод решения
        if (has_free) {
            std::cout << "System has infinitely many solutions:\n";
            int t_index = 1;
            for (int col = 0; col < count_variable; col++)
            {
                if (free_variables[col])
                {
                    std::cout << Variable[col] << " = t" << t_index++ << "\n";
                }
                else
                {
                    int r = -1;
                    for (int rr = 0; rr < row; rr++)
                    {
                        if (lead_col_of_row[rr] == col)
                        {
                            r = rr;
                            break;
                        }
                    }
                    std::cout << Variable[col] << " = " << free_coef[r];
                    for (int fc = 0; fc < count_variable; fc++)
                    {
                        if (free_variables[fc]) {
                            T coeff = -(*this)(r, fc);
                            if (std::abs(coeff) > epsilon) {
                                if (coeff > 0) {
                                    std::cout << " + ";
                                    if (std::abs(coeff - 1.0) > epsilon)
                                        std::cout << coeff;
                                }
                                else {
                                    std::cout << " - ";
                                    if (std::abs(coeff + 1.0) > epsilon)
                                        std::cout << -coeff;
                                }
                                std::cout << Variable[fc];
                            }
                        }
                    }
                    std::cout << "\n";
                }
            }
        }
        else {
            std::cout << "System solution:\n";
            for (int col = 0; col < count_variable; col++)
            {
                int r = -1;
                for (int rr = 0; rr < row; rr++)
                {
                    if (lead_col_of_row[rr] == col)
                    {
                        r = rr;
                        break;
                    }
                }
                std::cout << Variable[col] << " = " << free_coef[r] << "\n";
            }
        }
    }
};

int main()
{
    int count_variable;
    int count_strings;

    //количество переменных и их имена
    std::cout << "Enter the count variable\n";
    std::cin >> count_variable;
    vector<std::string> Variable(count_variable);
    for (int i = 0; i < count_variable; i++) {
        std::cout << "Enter the " << i + 1 << " variable name\n";
        std::cin >> Variable[i];
    }

    //количество строк и их ввод
    std::cout << "Enter the count strings\n";
    std::cin >> count_strings;
    Matrix<double> matrix(count_strings, count_variable);
    vector<double> free_coefficient(count_strings);
    for (int i = 0; i < count_strings; i++) {
        std::cout << "\nEnter the coefficients of variables the " << i + 1 << " strings" << std::endl;
        matrix.input_string(i);
        std::cout << "\nEnter the free coefficient of the " << i + 1 << " string\n";
        std::cin >> free_coefficient[i];
        std::cout << "\n";
    }

    matrix.print_matrix();
    try
    {
        matrix.solve_Gauss(free_coefficient, Variable);
    }
    catch (const std::exception& e)
    {
        std::cout << e.what() << std::endl;
    }
};