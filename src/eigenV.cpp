#include "eigenV.hpp"
#include "utility.hpp"
#include <vector>

using namespace arma;
using namespace corner;

bool eigenV::operator==(const eigenV& eV) const
{
    return (approx_equal(eigenvalue, eV.eigenvalue) && approx_equal(eigenvector, eV.eigenvector));
}

// Group eigenV with same value. Assuming that it is already sorted by abs value
void corner::groupWithSameAbs(std::vector<eigenV>& eV)
{
    check(eV.size() != 0, "groupEigenV", "passed empty vector");
    
    int i = 0;

    // Map containing the initial 
    while (i < eV.size())
    {
        double current_abs_eigenval = std::abs(eV[i].eigenvalue);

        // Count how many eigenvalues of equal abs are there
        int count_equal_abs = 0;
        for (int j = i; j != eV.size() && approx_equal(current_abs_eigenval, std::abs(eV[j].eigenvalue)); j++)
        {
            count_equal_abs++;
        }      

        int next_i = i + count_equal_abs;

        std::cout << "With abs: " << current_abs_eigenval << ": " << count_equal_abs << std::endl; 

        // Sort by real part
        std::sort(eV.begin() + i, eV.begin() + i + count_equal_abs, corner::sortByReal());

        groupWithSameReal(eV.begin() + i, eV.begin() + next_i);

        i = next_i;
    }

}

void corner::groupWithSameReal(std::vector<eigenV>::iterator beg, std::vector<eigenV>::iterator end)
{
  std::vector<eigenV>::iterator itr = beg;

    while (itr != end)
    {
        double current_real = std::real(itr->eigenvalue);

        // Count how many eigenvalues of equal real part are there
        int count_equal_real = 0;
        for (std::vector<eigenV>::iterator itr2 = itr; itr2 != end && approx_equal(current_real, std::real(itr2->eigenvalue)); itr2++)
        {
            count_equal_real++;
        }     

        std::sort(itr, itr + count_equal_real, sortByImag());

        itr = itr + count_equal_real;
    }

}


