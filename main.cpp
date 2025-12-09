
#include "FemObject.h"
#include <iostream>
#include <fstream>


int main(){




    // circleMesh(int n, int m, double r)
    // squareMesh(int n, double d, Eigen::Vector2d op)


    // BEGIN LAPLACE



    std::vector<int> n_levels_circle = {4, 6, 8, 10, 12, 14};
    // std::vector<int> m_levels_circle = {6, 8, 10, 12, 14, 16};

    int m = 6;
    std::vector<int> n_levels_square = {10, 14, 18, 22, 26, 30};


    std::ofstream laplaceFile("data/laplace/laplace.dat");
    std::ofstream poissonFile("data/poisson/poisson.dat");
    std::ofstream helmholtzFile("data/helmholtz/helmholtz.dat");

    laplaceFile << "#DoF\tError\n";
    poissonFile << "#DoF\tError\n";
    helmholtzFile << "#DoF\tError\n";


    FemObject fem;


    for(size_t k = 0; k < n_levels_circle.size(); k++){

        int n = n_levels_circle[k];



        fem.circleMesh(n,m,1.0);
        fem.setProblemType(LAPLACE);
        fem.solve();


        size_t DoF = fem.triang.getNodes()->size();
        double error = fem.getError();

        laplaceFile << DoF << "\t" << error << "\n";

        fem.visualization("data/laplace/laplace_n" + std::to_string(n) +
                          "_m" + std::to_string(m) + ".obj");

        fem.circleMesh(n, m, 1.0);
        fem.setProblemType(POISSON);
        fem.solve();

        DoF = fem.triang.getNodes()->size();
        error = fem.getError();

        poissonFile << DoF << "\t" << error << "\n";

        fem.visualization("data/poisson/poisson_n" + std::to_string(n) +
                          "_m" + std::to_string(m) + ".obj");

        // fem.squareMesh(40,1.0,{0.0,0.0});

    }


    for(int n: n_levels_square){
        fem.squareMesh(n, 1.0, {0.0,0.0});
        fem.setProblemType(HELMHOLTZ);
        fem.solve();

        auto DoF = fem.triang.getNodes()->size();
        auto error = fem.getError();

        helmholtzFile << DoF << "\t" << error << "\n";

        fem.visualization("data/helmholtz/helmholtz_n" + std::to_string(n) + ".obj");

    }




    // fem.circleMesh(8,10,1.0);

    // fem.setProblemType(LAPLACE);
    // fem.solve();

    // fem.visualization("laplace.obj");

    // auto DoF = fem.triang.getNodes()->size();



    // auto error = fem.getError();


    // WRITE OUT TO .DAT FILE SUCH THAT X IS DEGREES OF FREEMDOM #DoF and Y is L2 error, from variable error.



    // END LAPLACE





    // // BEGIN POISSON


    // fem.circleMesh(8,10,1.0);

    // fem.setProblemType(POISSON);

    // fem.solve();
    // fem.visualization("poisson.obj");


    // // END POISSON



    // // BEGIN HELMHOLTZ


    // fem.squareMesh(40,1.0,{0.0,0.0});
    // fem.setProblemType(HELMHOLTZ);
    // fem.solve();
    // fem.visualization("helmholtz.obj");



    // // END HELMHOLTZ


    // // BEGIN EIGENVALUE PROBLEM

    // fem.squareMesh(10,1.0,{0.0,0.0});
    // fem.setProblemType(EIGENVALUE);
    // fem.solve();
    // fem.visualization("eigenvalue.obj");

    // // END EIGENVALUE PROBLEM

    // auto error = fem.getError();
    // auto dof = fem.getDoF();

    return 0;
}
