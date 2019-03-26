#pragma once
#include <armadillo>

namespace corner {
    void check(const bool& condition, const char* class_name, const char* msg);
    void info(const bool& condition, const char* class_name, const char* msg);
    void warning(const bool& condition, const char* class_name, const char* msg);
    bool approx_equal(double, double, double tol=1E-8);
    bool approx_equal(arma::cx_double, arma::cx_double, double tol=1E-8);
    std::vector<std::string> split(const std::string& s, char delimiter);
    std::vector<std::string> splitString(std::string& s, std::string delimiter);
    class matrix
    {
        std::vector<std::vector<double> > m;
        public:
        matrix(const char*);

        public:
        double operator()(int i, int j=0);
        std::vector<double> col(int i);
        std::vector<double> row(int i);
    };
}