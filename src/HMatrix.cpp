#include "HMatrix.hpp"
#include "utility.hpp"

using namespace arma;
using namespace corner;

HMatrix::HMatrix(const cx_mat& m)
{
    check(m.n_rows == m.n_cols, "HMatrix::HMatrix", "Matrix must be square");
    this->m = m;
    
    // Get eigenvalues and eigenvectors in ascending order
    eig_gen(eigval, eigvec, m);
    uvec sorting_indices = sort_index(eigval);
    
    eigvec = eigvec.cols(sorting_indices);
    eigval = sort(eigval);
    
    //check(m*eigvec.col(size() - 1) == eigvec.col(size() - 1)*eigval(size() - 1), "HMatrix::HMatrix", "Eigval/Eigvec mismatch");
    //if(!approx_equal(m*eigvec.col(size() - 1), eigvec.col(size() - 1)*eigval(size() - 1), "absdiff", 1E-10))
   // {
     //   cout << m*eigvec.col(size() - 1) << endl << eigvec.col(size() - 1)*eigval(size() - 1) << endl;
    //}
    //check(approx_equal(m*eigvec.col(size() - 1), eigvec.col(size() - 1)*eigval(size() - 1), "absdiff", 1E-10), "HMatrix::HMatrix", "Eigval/Eigvec mismatch");
    
}

int HMatrix::size() const
{
    return m.n_rows;
}

arma::cx_mat HMatrix::GetONBasis()
{
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
            degvec.col(k - i) = eigvec.col(k);
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
    const cx_double lambda = eigval(i);
    int degeneration = 0;
    
    for (int k = 0; k < size(); k++)
    {
        if (corner::approx_equal(eigval(k),lambda))
            degeneration++;
    }
    //cout << "Degeneration: " << degeneration << endl;
    
    return degeneration;
}

cx_mat HMatrix::GramSchmidt(const cx_mat& deg_m)
{
    cx_mat ON_degvec(deg_m.n_rows, deg_m.n_cols);
    
    for (int k = 0; k < deg_m.n_cols; k++)
    {
        // Sum of projections of previous ON eigenvectors
        cx_vec proj_sum = zeros<cx_vec>(deg_m.n_rows);
        for (int p = 0; p < k ; p++)
        {
            proj_sum += Projection(deg_m.col(k), ON_degvec.col(p));
        }
        
        ON_degvec.col(k) = deg_m.col(k) - proj_sum;
    }
    
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
                warning(corner::approx_equal(cdot(ON_degvec.col(i), ON_degvec.col(j)),cx_double(1.,0.)), "HMatrix::GramSchmidt", "Vectors have no 1 norm");
                //cout << cdot(ON_degvec.col(i), ON_degvec.col(j)) << endl;
            }
            else
            {
                warning(corner::approx_equal(cdot(ON_degvec.col(i), ON_degvec.col(j)),cx_double(0.,0.)), "HMatrix::GramSchmidt", "Vectors are not normal");
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

const arma::cx_vec& HMatrix::GetEigenvalues() const
{
    return eigval;
}

bool HMatrix::IsHermitian()
{
    return arma::approx_equal(trans(m), m, "absdiff", 1E-10);
}

bool HMatrix::TraceOne()
{
    //cout << "Trace: " << trace(m) << endl;
    return corner::approx_equal(trace(m), cx_double(1.,0.));
}

bool HMatrix::EigSumIsOne()
{
    double sum=0.;
    for(int i=0; i<m.n_rows; i++)
    {
        sum+=real(eigval(i));
    }
    //cout << "Sum eigenvalues: " << sum << endl;
    return corner::approx_equal(sum, cx_double(1.,0.));
}

bool HMatrix::IsDM()
{
    return (HMatrix::IsHermitian() && HMatrix::TraceOne() && HMatrix::EigSumIsOne());
}

arma::cx_mat HMatrix::GetDM()
{
    cx_mat DM;
    int i=0;
    while(corner::approx_equal(eigval(i), cx_double(0.,0.)))
    {
        cx_vec rho_ss = eigvec.col(i);
        m = reshape(rho_ss, sqrt(size()), sqrt(size()));
        
        m = (m+trans(m))/2.;
        m = m/trace(m);
        HMatrix M(m);
         
        if(M.IsDM())
        {
            DM = m;
            break;
        }
        else
        {
            i+=1;
        }
    }
    
    return DM;
}