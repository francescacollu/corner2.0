/*
 Class describing an hermitian matrix
 */

#pragma once

#include<armadillo>

class HMatrix
{
    arma::cx_mat m;
    arma::cx_vec eigval;
    arma::cx_mat eigvec;
    bool IsLiouvillian;
    
    arma::cx_vec Projection(arma::cx_vec oldvec, arma::cx_vec newvec);

public:
    HMatrix(const arma::cx_mat&, bool IsLiouvillian=false);
    int size() const;
    
    // Build an orthonormal basis from the hermitian matrix
    arma::cx_mat GetONBasis();
    
    // Get order of degeneration of i-th eigenvalue
    int GetDegeneration(const int& i);
    
    // Apply Gram-Schmidt algorithm
    arma::cx_mat GramSchmidt(const arma::cx_mat&);
    
    //Get vector of sorted eigenvalues
    const arma::cx_vec& GetEigenvalues() const;
    
    // Is the matrix hermitian?
    bool IsHermitian() const;
    
    // Has the matrix trace one?
    bool TraceOne() const;
    
    // Has the matrix's eigenvalues sum equals to one?
    bool EigSumIsOne() const;

    // Has the matrix real and positive eigenvalues?
    bool EigRealPositive() const;
    
    // Has the matrix the properties of a density matrix?
    bool IsDM() const;

        // Is the matrix hermitian?
    static bool IsHermitian(const arma::cx_mat&);
    
    // Has the matrix trace one?
    static bool TraceOne(const arma::cx_mat&);
    
    // Has the matrix's eigenvalues sum equals to one?
    static bool EigSumIsOne(const arma::cx_mat&);

    // Has the matrix real and positive eigenvalues?
    static bool EigRealPositive(const arma::cx_mat&);
    
    // Has the matrix the properties of a density matrix?
    static bool IsDM(const arma::cx_mat&);

    // Calculates the steady state density matrix with the right parameters
    arma::cx_mat GetSteadyStateDM();
    
};
