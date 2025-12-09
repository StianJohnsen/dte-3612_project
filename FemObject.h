#ifndef FEMOBJECT_H
#define FEMOBJECT_H
#define _USE_MATH_DEFINES

#include <string>
#include <vector>
#include <ttl/halfedge/HeDart.h>
#include <ttl/halfedge/HeTraits.h>
#include <ttl/halfedge/HeTriang.h>

#include <Eigen/SparseCore>
#include <Eigen/Dense>
#include <Eigen/SparseCholesky>

#include <numeric>

#include <cmath>

enum ProblemType { LAPLACE, POISSON, HELMHOLTZ, EIGENVALUE };


class FemObject{

public:

    ProblemType problemType;
    hed::Triangulation triang;

    // Eigen::SparseMatrix<double> A; // Stiffness matrix
    // Eigen::SparseMatrix<double> M; // Mass matrix
    // Eigen::SparseMatrix<double> R; // Robin matrix
    // Eigen::VectorXd b; // Load vector
    // Eigen::VectorXd r; // Robin vector


    void setProblemType(ProblemType problemType);
    void solve();

    void circleMesh(int n, int m, double r);
    void squareMesh(int n, double d, Eigen::Vector2d Op);

    Eigen::SparseMatrix<double> massMat(const std::list<hed::Edge*>& leading_edges, int np,const std::unordered_map<const hed::Node*, int>& nodeToIndex);
    Eigen::SparseMatrix<double> stiffMat(const std::list<hed::Edge*>& leading_edges, int np, const std::unordered_map<const hed::Node*, int>& nodeToIndex);
    Eigen::VectorXd loadVect(const std::list<hed::Edge*>& leading_edges, int np, const std::unordered_map<const hed::Node*, int>& nodeToIndex);
    Eigen::SparseMatrix<double> robinMat(const std::list<hed::Dart>& boundary, int np, const std::unordered_map<const hed::Node*, int>& nodeToIndex);
    Eigen::VectorXd robinVect(const std::list<hed::Dart>& boundary, int np, const std::unordered_map<const hed::Node*, int>& nodeToIndex);

    double kappa(double x, double y);

    // double kappa(double x, double y); // BC type
    double gN(double x, double y); // Neumann BC
    double gD(double x, double y); // Dirichlet BC
    double f(double x, double y); // SOURCE FUNCTION
    double u_exact(double x, double y); // Exact solution

    double triArea(hed::Node* N1,hed::Node* N2,hed::Node* N3);
    Eigen::Vector2<Eigen::Vector3d> gradients(Eigen::Vector3d x, Eigen::Vector3d y, double area);

    // double getError(const std::list<hed::Edge*>& leading_edges,const std::unordered_map<const hed::Node*, int>& nodeToIndex); // Compare approximated solution to the exact solution

    // double getError(const std::list<hed::Edge*>& leading_edges); // Compare approximated solution to the exact solution

    double getError();

    int getDoF(); // Number of nodes // Degrees of Freedom

    void visualization(std::string filename); // Create .obj file

    std::string problemTypeToString(ProblemType pt);






};

#endif // FEMOBJECT_H
