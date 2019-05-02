#include "HMatrix.hpp"
#include "utility.hpp"
#include <algorithm>

using namespace arma;
using namespace corner;

HMatrix::HMatrix(const cx_mat& m, bool IsLiouvillian)
{
    check(m.n_rows == m.n_cols, "HMatrix::HMatrix", "Matrix must be square");
    this->m = m;
    this->IsLiouvillian = IsLiouvillian;
    cx_vec eigval;
    cx_mat eigvec;

    // Get eigenvalues and eigenvectors in ascending order
    eig_gen(eigval, eigvec, m);

    for(int i = 0; i < eV.size(); i++)
    {
       check(approx_equal(m*(eigvec.col(i)), eigvec.col(i)*eigval(i), "absdiff", 1E-3), "HMatrix::HMatrix", "Pre-sort Eigval/Eigvec mismatch");
       //cout << i << endl;
    }

    for(int i = 0; i < eigval.n_rows; i++)
    {
        corner::eigenV e;
        e.eigenvalue = eigval(i);
        e.eigenvector = eigvec.col(i);
        eV.push_back(e);
    }

    // Sort in lambda absolute value 
    // Problem: it could sort it like this: 1 -1 1 1 2 ...    
    std::sort(eV.begin(), eV.end(), sortByAbs());

    // Group eigenvalues with same value to solve the problem
    corner::groupWithSameAbs(eV);



    for(int i = 0; i < eV.size(); i++)
    {
        check(approx_equal(m*(eV[i].eigenvector), eV[i].eigenvector*eV[i].eigenvalue, "absdiff", 1E-3), "HMatrix::HMatrix", "Eigval/Eigvec mismatch");
        //cout << i << endl;
    }

    int i = 0;
    while (i < size())
    {
        int degeneration = GetDegeneration(i);

        cx_double lambda = eV[i].eigenvalue;

        for (int j = 0; j < degeneration; j++)
        {
            check(approx_equal(eV[i+j].eigenvalue, lambda), "HMatrix::HMatrix", "Degenerate eigvals are not aligned");
            check(approx_equal(m*eV[i+j].eigenvector, eV[i+j].eigenvector*lambda, "absdiff", 1E-3), "HMatrix::HMatrix", "Eigval/Eigvec mismatch");
        }

        i += degeneration;
    }

    if(!IsLiouvillian)
    {
        check(EigRealPositive(), "HMatrix::HMatrix", "eigenvalues not real or positive");
    }
}

int HMatrix::size() const
{
    return m.n_rows;
}

arma::cx_mat HMatrix::GetONBasis()
{
    check(!IsLiouvillian, "HMatrix::GetONBasis", "Cannot call this method for a liouvillian matrix");
    // Basis matrix
    cx_mat B(size(), size());
    
    int i = 0;
    while (i < size())
    {
        int degeneration = GetDegeneration(i);
        
        // Matrix of eigenvectors corresponding to the same eigenvalue
        cx_mat degvec(size(), degeneration);
        
        for (int k = i; k < i + degeneration; k++)
        {
            degvec.col(k - i) = eV[k].eigenvector;
        }
        
        // Orthonormalization via Gram-Schmidt algorithm
        cx_mat ON_degvec = GramSchmidt(degvec);
        
        B.submat(0, i, size()-1, i + degeneration-1) = ON_degvec;
        
        // Get ready for next eigenvalue
        i += degeneration;
    }
    
    // Assert if basis is ON
    for (int i = 0; i < B.n_cols; i++)
    {
        for (int j = 0; j < B.n_cols; j++)
        {
            if (i == j)
            {
                warning(corner::approx_equal(cdot(B.col(i), B.col(j)), cx_double(1., 0.)), "HMatrix::GetONBasis", "Vectors have no 1 norm");
                //cout << cdot(B.col(i), B.col(j)) << endl;
            }
            else
            {
                warning(corner::approx_equal(cdot(B.col(i), B.col(j)), cx_double(0.,0.)), "HMatrix::GetONBasis", "Vectors are not normal");
                if(!corner::approx_equal(cdot(B.col(i), B.col(j)), cx_double(0.,0.)))
                {
                    cout << cdot(B.col(i), B.col(j)) << endl;
                }
            }
        }
    }
    
    return B;
}

int HMatrix::GetDegeneration(const int& i)
{
    const cx_double lambda = eV[i].eigenvalue;
    int degeneration = 0;
    
    for (int k = 0; k < size(); k++)
    {
        if (corner::approx_equal(eV[k].eigenvalue, lambda))
            degeneration++;
    }
    
    return degeneration;
}

cx_mat HMatrix::GramSchmidt(const cx_mat& deg_m)
{
    check(rank(deg_m) == std::min(deg_m.n_rows, deg_m.n_cols), "HMatrix::GramSchmidt", "Vectors not l.i.");
    cx_mat ON_degvec(deg_m.n_rows, deg_m.n_cols);

    for (int k = 0; k < deg_m.n_cols; k++)
    {
        ON_degvec.col(k) = deg_m.col(k);
        for (int p = 0; p < k ; p++)
        {
            ON_degvec.col(k) = ON_degvec.col(k) - Projection(ON_degvec.col(k), ON_degvec.col(p));
        }
        ON_degvec.col(k) = ON_degvec.col(k) / sqrt(cdot(ON_degvec.col(k),ON_degvec.col(k)));
    }
    
    // for (int k = 0; k < deg_m.n_cols; k++)
    // {
    //     // Sum of projections of previous ON eigenvectors
    //     cx_vec proj_sum = zeros<cx_vec>(deg_m.n_rows);
    //     for (int p = 0; p < k ; p++)
    //     {
    //         proj_sum += Projection(deg_m.col(k), ON_degvec.col(p));
    //     }
        
    //     ON_degvec.col(k) = deg_m.col(k) - proj_sum;
    // }
    
    check(ON_degvec.n_cols == deg_m.n_cols, "HMatrix::GramSchmidt", "New vector size is wrong");
    
    // Normalize basis
    for (int k = 0; k < ON_degvec.n_cols; k++)
    {
        //cout << "norm(ON_degvec(first): " <<  sqrt(cdot(ON_degvec.col(k),ON_degvec.col(k))) << endl;
        ON_degvec.col(k) = ON_degvec.col(k) / sqrt(cdot(ON_degvec.col(k),ON_degvec.col(k)));
        //cout << "norm(ON_degvec(then): " <<  sqrt(cdot(ON_degvec.col(k),ON_degvec.col(k))) << endl;
        
    }
    
    // Assert if matrix is ON
    for (int i = 0; i < ON_degvec.n_cols; i++)
    {
        for (int j = 0; j < ON_degvec.n_cols; j++)
        {
            if (i == j)
            {
                check(corner::approx_equal(cdot(ON_degvec.col(i), ON_degvec.col(j)),cx_double(1.,0.)), "HMatrix::GramSchmidt", "Vectors have no 1 norm");
                //cout << cdot(ON_degvec.col(i), ON_degvec.col(j)) << endl;
            }
            else
            {
                check(corner::approx_equal(cdot(ON_degvec.col(i), ON_degvec.col(j)),cx_double(0.,0.)), "HMatrix::GramSchmidt", "Vectors are not normal");
                //cout << cdot(ON_degvec.col(i), ON_degvec.col(j)) <<endl;
            }
        }
    }
    return ON_degvec;
}

arma::cx_vec HMatrix::Projection(cx_vec oldvec, cx_vec newvec)
{
    return (cdot(newvec, oldvec) / cdot(newvec,newvec))*newvec;
}

arma::cx_vec HMatrix::GetEigenvalues() const
{
    cx_vec v;
    for (int i = 0; i < eV.size(); i++)
    {
        v(i) = eV[i].eigenvalue;
    }
    return v;
}

bool HMatrix::IsHermitian() const
{
    return arma::approx_equal(trans(m), m, "absdiff", 1E-10);
}

bool HMatrix::TraceOne() const
{
    //cout << "Trace: " << trace(m) << endl;
    return corner::approx_equal(trace(m), cx_double(1.,0.));
}

bool HMatrix::EigSumIsOne() const
{
    double sum=0.;
    for(int i=0; i<m.n_rows; i++)
    {
        sum+=real(eV[i].eigenvalue);
    }
    //cout << "Sum eigenvalues: " << sum << endl;
    return corner::approx_equal(sum, cx_double(1.,0.));
}

bool HMatrix::EigRealPositive() const
{
    for(int i=0; i<m.n_rows; i++)
    {
        if(!approx_equal(imag(eV[i].eigenvalue), 0.0, 1E-10) || real(eV[i].eigenvalue) < -1E-10)
        {
            std::cout << "Found unsuitable eigenvalue: " << eV[i].eigenvalue << std::endl;
            return false;
        }
    }
    return true;
}

bool HMatrix::IsDM() const
{
    return (IsHermitian() && TraceOne() && EigSumIsOne() && EigRealPositive());
}

bool HMatrix::IsHermitian(const cx_mat& m)
{
    return arma::approx_equal(trans(m), m, "absdiff", 1E-10);
}

bool HMatrix::TraceOne(const cx_mat& m)
{
    //cout << "Trace: " << trace(m) << endl;
    return corner::approx_equal(trace(m), cx_double(1.,0.));
}

bool HMatrix::EigSumIsOne(const cx_mat& m)
{
    double sum=0.;
    cx_vec eigval;
    eig_gen(eigval, m);
    for(int i=0; i<m.n_rows; i++)
    {
        sum+=real(eigval(i));
    }
    //cout << "Sum eigenvalues: " << sum << endl;
    return corner::approx_equal(sum, cx_double(1.,0.));
}

bool HMatrix::EigRealPositive(const cx_mat& m) 
{
    cx_vec eigval;
    eig_gen(eigval, m);
    for(int i=0; i<m.n_rows; i++)
    {
        if(!approx_equal(imag(eigval(i)), cx_double(0., 0.), 1E-10) || real(eigval(i)) < -1E-10)
        {
            return false;
        }
    }
    return true;
}

bool HMatrix::IsDM(const cx_mat& m)
{
    return (IsHermitian(m) && TraceOne(m) && EigSumIsOne(m) && EigRealPositive(m));
}

arma::cx_mat HMatrix::GetSteadyStateDM()
{
    //std::vector<cx_vec> eigvecs;
    cx_mat dm;

    //cout << "Degeneration: " << GetDegeneration(0) << endl;

    check(corner::approx_equal(eV[0].eigenvalue, arma::cx_double(0,0)), "HMatrix::GetSteadyStateDM", "Null eigenvalue not found.");

    // std::vector<cx_mat> eigvecs;
    // for(int i=0; i<GetDegeneration(0); i++)
    // {
    //    eigvecs.push_back(eV[i].eigenvector);
    // }

    // cx_vec eigvecTot = 0.*eigvecs[0];
    // for(int i=0; i<GetDegeneration(0); i++)
    // {
    //    eigvecTot += (1/float(GetDegeneration(0)))*eigvecs[i];
    // }

    bool found = false;
    for (int i = 0; i < GetDegeneration(0); i++)
    {
        dm = reshape(eV[i].eigenvector, sqrt(size()), sqrt(size()));
        dm = (dm+trans(dm))/2.;
        dm = dm/trace(dm);

        if (HMatrix::IsDM(dm))
        {
            found = true;
            break;
        }
    }
    
    check(found, "HMatrix::GetSteadyStateDM", "Did not find a suitable DM");

    check(!corner::approx_equal(cx_double(trace(dm)), cx_double(0.,0.)), "HMatrix::GetSteadyStateDM", "trace(dm) approximately equals to zero.");
    
    check(HMatrix::IsDM(dm), "HMatrix::GetSteadyStateDM", "The matrix does not satisfy one ore more than one property of the DM. A check is necessessary.");

    return dm;
}