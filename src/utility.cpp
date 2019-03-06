#include "utility.hpp"
#include <iostream>
#include <exception>

void corner::check(const bool& condition, const char* class_name, const char* msg)
{
    if (!condition)
    {
        std::cout << "ERROR in " << class_name << " : " << msg << std::endl;
        throw std::exception();
    }
}

void corner::info(const bool& condition, const char* class_name, const char* msg)
{
    if (!condition)
    {
        std::cout << "INFO in " << class_name << " : " << msg << std::endl;
    }
}

void corner::warning(const bool& condition, const char* class_name, const char* msg)
{
    if (!condition)
    {
        std::cout << "WARNING in " << class_name << " : " << msg << std::endl;
    }
}

bool corner::approx_equal(double val1, double val2, double tol )
{
    //std::cout << "double_approx\n";
    if(std::abs(val1-val2)<tol)
    {
        return true;
    }
    return false;
}

bool corner::approx_equal(arma::cx_double val1, arma::cx_double val2, double tol )
{
    //std::cout << "cx_approx" << std::endl;
    if(std::abs(std::sqrt(conj(val1-val2)*(val1-val2)))<tol)
    {
        return true;
    }
    return false;
}
