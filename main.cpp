
#include "FemObject.h"
#include <iostream>


int main(){

    FemObject fem;




    // // BEGIN LAPLACE



    // fem.circleMesh(8,10,1.0);
    // // fem.squareMesh(10,10.0,{0.0,0.0});

    // fem.setProblemType(LAPLACE);
    // fem.solve();

    // fem.visualization("laplace.obj");



    // // END LAPLACE





    // // BEGIN POISSON


    // fem.circleMesh(8,10,1.0);

    // fem.setProblemType(POISSON);

    // fem.solve();
    // fem.visualization("poisson.obj");


    // // END POISSON



    // BEGIN HELMHOLTZ


    // fem.squareMesh(10,10.0,{0.0,0.0});
    // fem.setProblemType(HELMHOLTZ);
    // fem.solve();
    // fem.visualization("helmholtz.obj");



    // END HELMHOLTZ


    // BEGIN EIGENVALUE PROBLEM

    fem.squareMesh(10,10.0,{0.0,0.0});
    fem.setProblemType(EIGENVALUE);
    fem.solve();
    fem.visualization("eigenvalue.obj");

    // END EIGENVALUE PROBLEM

    // auto error = fem.getError();
    // auto dof = fem.getDoF();

    return 0;
}
