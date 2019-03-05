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
    
    arma::cx_vec Projection(arma::cx_vec oldvec, arma::cx_vec newvec);

public:
    HMatrix(const arma::cx_mat&);
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

    // Calculates the steady state density matrix with the right parameters
    arma::cx_mat GetSteadyStateDM();
    
};
