#pragma once

#include <exception>
#include <string>
#include <sstream>
#include "vector.h"

class IncompatibleVectorException : std::exception
{
    size_t size1 = 0, size2 = 0;
    public:
    IncompatibleVectorException(const corner::vector& v1, const corner::vector& v2)
    {
        size1 = v1.size();
        size2 = v2.size();
    }

    const char* what()
    {
        std::stringstream ss;
        std::string s1, s2;
        ss << size1;
        ss >> s1;
        ss.clear();

        ss << size2;
        ss >> s2;

        return ((std::string("Incompatible vector sizes: ") + s1 + " and " + s2).c_str());
    }
};