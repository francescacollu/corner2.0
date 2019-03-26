#include "utility.hpp"
#include <iostream>
#include <exception>
#include <string>

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

std::vector<std::string> corner::split(const std::string& s, char delimiter)
{
   std::vector<std::string> tokens;
   std::string token;
   std::istringstream tokenStream(s);
   while (std::getline(tokenStream, token, delimiter))
   {
      tokens.push_back(token);
   }
   return tokens;
}

std::vector<std::string> corner::splitString(std::string& s, std::string delimiter)
{
    std::vector<std::string> tokens;
    std::string token;
    size_t pos = 0;
    while ((pos = s.find(delimiter)) != std::string::npos) 
    {
        token = s.substr(0, pos);
        tokens.push_back(token);
        s.erase(0, pos + delimiter.length());
    }

    return tokens;
}

double corner::matrix::operator()(int i, int j)
{
    corner::check(i >= 0 && j >= 0, "corner::matrix::operator()", "Not valid indices.\n");
    return m[i][j];
}

std::vector<double> corner::matrix::col(int j)
{
    corner::check(m.size()!=0, "corner::matrix::col()", "matrix is empty.");
    corner::check(j>=0 && j<m[0].size(), "corner::matrix::col()", "invalid index.");
    std::vector<double> vec;
    for(int i = 0; i < m.size(); i++)
    {
        vec.push_back(m[i][j]);
    }
    return vec;
}

std::vector<double> corner::matrix::row(int i)
{
    corner::check(i>=0 && i<m.size(), "corner::matrix::row()", "invalid index.");
    return m[i];
}

corner::matrix::matrix(const char* c)
{
    std::ifstream filename(c);
    std::string s;

    corner::check(filename.is_open(), "corner::matrix::matrix", "The file has not been opened");

    while(!filename.eof()){
        std::getline(filename, s);
        std::vector<std::string> vec = corner::split(s, '\t');
        if(vec.size() == 0) continue;
        if(m.size()!=0)
        {
            corner::check(m.back().size() == vec.size(), "corner::matrix::matrix", "There is a problem with the matrix dimension.");
        }
        std::vector<double> newVec;
        for(int k = 0; k < vec.size(); k++)
        {   
            newVec.push_back(atof(vec[k].c_str()));
        }
        m.push_back(newVec);
    }
}