
#include "FemObject.h"
#include <iostream>


int main(){

    FemObject fem;
    // fem.circleMesh(8,10,1.0);
    fem.squareMesh(10,10.0,{0.0,0.0});

    // fem.setProblemType(EIGENVALUE);
    // fem.solve()

    // fem.visualization("square.obj");

    // auto error = fem.getError();
    // auto dof = fem.getDoF();

    return 0;
}
