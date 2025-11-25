#include "FemObject.h"


void FemObject::circleMesh(int n, int m, double r){

    std::vector<hed::Node*> nodes;

    nodes.reserve(1 + n * m * (n + 1) / 2);
    nodes.push_back(new hed::Node(0.0,0.0));

    for(size_t j = 0; j < n; ++j){

        double radius = (static_cast<double>(j+1) / n) * r;

        int numPoints = m * (j + 1);

        for(size_t i = 0; i < numPoints ; ++i){

            double alpha = (2.0 * M_PI * i) / numPoints;
            double x = radius * std::cos(alpha);
            double y = radius * std::sin(alpha);

            nodes.push_back(new hed::Node(x,y));

        }

    }

    this->triang.createDelaunay(nodes.begin(),nodes.end());
};

void FemObject::squareMesh(int n, double d, Eigen::Vector2d Op){

    std::vector<hed::Node*> nodes;

    nodes.reserve(n*n);


    for(size_t j = 0; j < n; ++j){

        double y = Op.y() + (static_cast<double>(j) / (n-1)) * d;

        for(size_t i = 0; i < n; ++i){

            double x = Op.x() + static_cast<double>(i) / (n-1) * d;

            nodes.push_back(new hed::Node(x,y));
        }

    }

    this->triang.createDelaunay(nodes.begin(),nodes.end());

};

void FemObject::visualization(std::string filename) {

    std::ofstream obj_file(filename);
    if (!obj_file.is_open()) {
        std::cerr << "Failed to open " << filename << "\n";
        return;
    }

    // --------------------------------------------------------------------
    // 1. Get all nodes from the triangulation
    // --------------------------------------------------------------------
    const std::list<hed::Node*>* nodes = this->triang.getNodes();  // returns pointer
    if (!nodes) {
        std::cerr << "No nodes found in triangulation.\n";
        return;
    }

    std::unordered_map<const hed::Node*, int> nodeIndex;
    int idx = 1;

    for (const auto& node : *nodes) {
        obj_file << "v " << node->x() << " " << node->y() << " " << node->z() << "\n";
        nodeIndex[node] = idx++;
    }

    // --------------------------------------------------------------------
    // 2. Iterate through all leading edges (one per triangle)
    // --------------------------------------------------------------------
    const std::list<hed::Edge*>& leading_edges = this->triang.getLeadingEdges();

    for (const auto& edge : leading_edges) {
        // Each leading edge represents one triangle (CCW orientation)
        hed::Edge* e1 = edge;
        hed::Edge* e2 = e1->getNextEdgeInFace();
        hed::Edge* e3 = e2->getNextEdgeInFace();

        if (!e1 || !e2 || !e3) continue; // safety check

        hed::Node* n1 = e1->getSourceNode();
        hed::Node* n2 = e2->getSourceNode();
        hed::Node* n3 = e3->getSourceNode();

        if (!n1 || !n2 || !n3) continue; // skip malformed triangles

        // Write a face (indices are 1-based in OBJ)
        obj_file << "f "
                 << nodeIndex[n1] << " "
                 << nodeIndex[n2] << " "
                 << nodeIndex[n3] << "\n";
    }

    obj_file.close();
    std::cout << "Mesh exported to " << filename <<  " (" << nodes->size()
              << " vertices, " << leading_edges.size() << " triangles)\n";
}





double FemObject::triArea(hed::Node* N1,hed::Node* N2,hed::Node* N3){

    Eigen::Vector3d a = {N2->x() - N1->x() , N2->y() - N1->y() ,0};
    Eigen::Vector3d b = {N3->x() - N1->x(),N3->y() - N1->y(),0};
    double area = 0.5 * (b.cross(a)).norm();
    return area;

}


Eigen::Vector2<Eigen::Vector3d> FemObject::gradients(Eigen::Vector3d x, Eigen::Vector3d y, double area){

    Eigen::Vector3d b;
    b(0) = (y(1)-y(2)) / (2*area);
    b(1) = (y(2)-y(0)) / (2*area);
    b(2) = (y(0)-y(1)) / (2*area);

    Eigen::Vector3d c;
    c(0) = (x(2)-x(1)) / (2*area);
    c(1) = (x(0)-x(2)) / (2*area);
    c(2) = (x(1)-x(0)) / (2*area);

    return {b, c};

}



Eigen::SparseMatrix<double> FemObject::stiffMat(const std::list<hed::Edge*>& leading_edges, int np){

    Eigen::SparseMatrix<double> A;
    A.resize(np,np);
    A.reserve(Eigen::VectorXi::Constant(np,10));


    // const std::list<hed::Edge*>& leading_edges = this->triang.getLeadingEdges();

    for(const auto& edge : leading_edges){
        hed::Edge* e1 = edge;
        hed::Edge* e2 = e1->getNextEdgeInFace();
        hed::Edge* e3 = e2->getNextEdgeInFace();

        if (!e1 || !e2 || !e3) continue; // safety check

        hed::Node* n1 = e1->getSourceNode();
        hed::Node* n2 = e2->getSourceNode();
        hed::Node* n3 = e3->getSourceNode();

        if (!n1 || !n2 || !n3) continue; // skip malformed triangles

        Eigen::Vector3i loc2glb = {n1->id(),n2->id(),n3->id()};


        Eigen::Vector3d x_coords = {n1->x(),n2->x(),n3->x()};
        Eigen::Vector3d y_coords = {n1->y(),n2->y(),n3->y()};
        double area = triArea(n1,n2,n3);

        Eigen::Vector2<Eigen::Vector3d> gradPhi = gradients(x_coords,y_coords,area);

        Eigen::Vector3d b = gradPhi(0);
        Eigen::Vector3d c = gradPhi(1);

        Eigen::Matrix3d AK = (b*b.transpose() + c*c.transpose()) * area;

        for(int i = 0; i < 3; ++i){

            int I = loc2glb(i);

            for(int j = 0; j < 3; ++j){

                int J = loc2glb(j);

                A.coeffRef(I,J) += AK(i,j);

            }
        }

    }

    A.makeCompressed();
    return A;

}



Eigen::SparseMatrix<double> FemObject::massMat(const std::list<hed::Edge*>& leading_edges, int np){

    Eigen::SparseMatrix<double> M;
    M.resize(np,np);
    M.reserve(Eigen::VectorXi::Constant(np,10));

    for(const auto& edge : leading_edges){
        hed::Edge* e1 = edge;
        hed::Edge* e2 = e1->getNextEdgeInFace();
        hed::Edge* e3 = e2->getNextEdgeInFace();

        if (!e1 || !e2 || !e3) continue; // safety check

        hed::Node* n1 = e1->getSourceNode();
        hed::Node* n2 = e2->getSourceNode();
        hed::Node* n3 = e3->getSourceNode();

        if (!n1 || !n2 || !n3) continue; // skip malformed triangles

        Eigen::Vector3i loc2glb = {n1->id(),n2->id(),n3->id()};


        double area = triArea(n1,n2,n3);


        Eigen::Matrix3d MK;
        MK <<  2,1,1,
                1,2,1,
                1,1,2;

        MK *= (area / 12.0);

        for(int i = 0; i < 3; ++i){

            int I = loc2glb(i);

            for(int j = 0; j < 3; ++j){

                int J = loc2glb(j);

                M.coeffRef(I,J) += MK(i,j);
            }
        }
    }

    M.makeCompressed();
    return M;

}


Eigen::VectorXd FemObject::loadVect(const std::list<hed::Edge*>& leading_edges, int np){

    Eigen::VectorXd b = Eigen::VectorXd::Zero(np);


    for(const auto& edge : leading_edges){
        hed::Edge* e1 = edge;
        hed::Edge* e2 = e1->getNextEdgeInFace();
        hed::Edge* e3 = e2->getNextEdgeInFace();

        if (!e1 || !e2 || !e3) continue; // safety check

        hed::Node* n1 = e1->getSourceNode();
        hed::Node* n2 = e2->getSourceNode();
        hed::Node* n3 = e3->getSourceNode();

        if (!n1 || !n2 || !n3) continue; // skip malformed triangles

        Eigen::Vector3i loc2glb = {n1->id(),n2->id(),n3->id()};


        Eigen::Vector3d x_coords = {n1->x(),n2->x(),n3->x()};
        Eigen::Vector3d y_coords = {n1->y(),n2->y(),n3->y()};
        double area = triArea(n1,n2,n3);

        double x_c = (x_coords(0)+x_coords(1)+x_coords(2)) / 3;
        double y_c = (y_coords(0)+y_coords(1)+y_coords(2)) / 3;

        double F = f(x_c,y_c, ProblemType::EIGEN);

        Eigen::Vector3d vec = {1,1,1};
        Eigen::Vector3d bK = F / 3 * vec * area;

        for(int j = 0; j < 3; ++j){

            int J = loc2glb(j);

            b.coeffRef(J) += bK(j);
        }


    }

    return b;
}


double FemObject::f(double x, double y, ProblemType problemType, double k){

    switch (problemType ) {
    case LAPLACE:
        return 0.0;

    case POISSON:

        return (2.0 * M_PI * M_PI + k*k)
               * std::sin(M_PI * x) * std::sin(M_PI * y);

    case EIGEN:
        return 0.0;

    default:
        return 0.0;
    }

}

