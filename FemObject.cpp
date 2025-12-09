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

        float y = Op.y() + (static_cast<double>(j) / (n-1)) * d;

        for(size_t i = 0; i < n; ++i){

            float x = Op.x() + static_cast<double>(i) / (n-1) * d;

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
        obj_file << "v " << node->x() << " " << node->z() << " " << node->y() << "\n";
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
                 << nodeIndex[n3] << " "
                 << nodeIndex[n2] << "\n";
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



Eigen::SparseMatrix<double> FemObject::stiffMat(
    const std::list<hed::Edge*>& leading_edges,
    int np,
    const std::unordered_map<const hed::Node*, int>& nodeToIndex)
{
    Eigen::SparseMatrix<double> A(np, np);
    A.reserve(Eigen::VectorXi::Constant(np, 6));

    for (auto edge : leading_edges) {
        hed::Edge* e1 = edge;
        hed::Edge* e2 = e1->getNextEdgeInFace();
        hed::Edge* e3 = e2->getNextEdgeInFace();

        hed::Node* n1 = e1->getSourceNode();
        hed::Node* n2 = e2->getSourceNode();
        hed::Node* n3 = e3->getSourceNode();

        int I1 = nodeToIndex.at(n1);
        int I2 = nodeToIndex.at(n2);
        int I3 = nodeToIndex.at(n3);

        Eigen::Vector3d x = {n1->x(), n2->x(), n3->x()};
        Eigen::Vector3d y = {n1->y(), n2->y(), n3->y()};

        double area = triArea(n1, n2, n3);
        auto grad = gradients(x, y, area);

        Eigen::Vector3d b = grad(0);
        Eigen::Vector3d c = grad(1);

        Eigen::Matrix3d AK = (b * b.transpose() + c * c.transpose()) * area;

        int idx[3] = {I1, I2, I3};

        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                A.coeffRef(idx[i], idx[j]) += AK(i, j);
    }

    return A;
}


// Eigen::SparseMatrix<double> FemObject::stiffMat(const std::list<hed::Edge*>& leading_edges, int np){

//     Eigen::SparseMatrix<double> A;
//     A.resize(np,np);
//     A.reserve(Eigen::VectorXi::Constant(np,10));

//     // const std::list<hed::Edge*>& leading_edges = this->triang.getLeadingEdges();

//     for(const auto& edge : leading_edges){
//         hed::Edge* e1 = edge;
//         hed::Edge* e2 = e1->getNextEdgeInFace();
//         hed::Edge* e3 = e2->getNextEdgeInFace();

//         if (!e1 || !e2 || !e3) continue; // safety check

//         hed::Node* n1 = e1->getSourceNode();
//         hed::Node* n2 = e2->getSourceNode();
//         hed::Node* n3 = e3->getSourceNode();

//         if (!n1 || !n2 || !n3) continue; // skip malformed triangles

//         Eigen::Vector3i loc2glb = {n1->id(),n2->id(),n3->id()};


//         Eigen::Vector3d x_coords = {n1->x(),n2->x(),n3->x()};
//         Eigen::Vector3d y_coords = {n1->y(),n2->y(),n3->y()};
//         double area = triArea(n1,n2,n3);

//         Eigen::Vector2<Eigen::Vector3d> gradPhi = gradients(x_coords,y_coords,area);

//         Eigen::Vector3d b = gradPhi(0);
//         Eigen::Vector3d c = gradPhi(1);

//         Eigen::Matrix3d AK = (b*b.transpose() + c*c.transpose()) * area;

//         for(int i = 0; i < 3; ++i){

//             int I = loc2glb(i);

//             for(int j = 0; j < 3; ++j){

//                 int J = loc2glb(j);

//                 A.coeffRef(I,J) += AK(i,j);

//             }
//         }

//     }

//     // A.makeCompressed();
//     return A;

// }


Eigen::SparseMatrix<double> FemObject::massMat(
    const std::list<hed::Edge*>& leading_edges,
    int np,
    const std::unordered_map<const hed::Node*, int>& nodeToIndex)
{
    Eigen::SparseMatrix<double> M(np, np);
    M.reserve(Eigen::VectorXi::Constant(np, 6));

    for (auto edge : leading_edges) {
        hed::Edge* e1 = edge;
        hed::Edge* e2 = e1->getNextEdgeInFace();
        hed::Edge* e3 = e2->getNextEdgeInFace();

        hed::Node* n1 = e1->getSourceNode();
        hed::Node* n2 = e2->getSourceNode();
        hed::Node* n3 = e3->getSourceNode();

        int I1 = nodeToIndex.at(n1);
        int I2 = nodeToIndex.at(n2);
        int I3 = nodeToIndex.at(n3);

        double area = triArea(n1, n2, n3);

        Eigen::Matrix3d MK;
        MK << 2, 1, 1,
            1, 2, 1,
            1, 1, 2;
        MK *= area / 12.0;

        int idx[3] = {I1, I2, I3};

        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                M.coeffRef(idx[i], idx[j]) += MK(i, j);
    }

    return M;
}



// Eigen::SparseMatrix<double> FemObject::massMat(const std::list<hed::Edge*>& leading_edges, int np){

//     Eigen::SparseMatrix<double> M;
//     M.resize(np,np);
//     M.reserve(Eigen::VectorXi::Constant(np,10));

//     for(const auto& edge : leading_edges){
//         hed::Edge* e1 = edge;
//         hed::Edge* e2 = e1->getNextEdgeInFace();
//         hed::Edge* e3 = e2->getNextEdgeInFace();

//         if (!e1 || !e2 || !e3) continue; // safety check

//         hed::Node* n1 = e1->getSourceNode();
//         hed::Node* n2 = e2->getSourceNode();
//         hed::Node* n3 = e3->getSourceNode();

//         if (!n1 || !n2 || !n3) continue; // skip malformed triangles

//         Eigen::Vector3i loc2glb = {n1->id(),n2->id(),n3->id()};


//         double area = triArea(n1,n2,n3);


//         Eigen::Matrix3d MK;
//         MK <<  2,1,1,
//                 1,2,1,
//                 1,1,2;

//         MK *= (area / 12.0);

//         for(int i = 0; i < 3; ++i){

//             int I = loc2glb(i);

//             for(int j = 0; j < 3; ++j){

//                 int J = loc2glb(j);

//                 M.coeffRef(I,J) += MK(i,j);
//             }
//         }
//     }

//     // M.makeCompressed();
//     return M;

// }

Eigen::VectorXd FemObject::loadVect(
    const std::list<hed::Edge*>& leading_edges,
    int np,
    const std::unordered_map<const hed::Node*, int>& nodeToIndex)
{
    Eigen::VectorXd b = Eigen::VectorXd::Zero(np);

    for (auto edge : leading_edges) {
        hed::Edge* e1 = edge;
        hed::Edge* e2 = e1->getNextEdgeInFace();
        hed::Edge* e3 = e2->getNextEdgeInFace();

        hed::Node* n1 = e1->getSourceNode();
        hed::Node* n2 = e2->getSourceNode();
        hed::Node* n3 = e3->getSourceNode();

        int I1 = nodeToIndex.at(n1);
        int I2 = nodeToIndex.at(n2);
        int I3 = nodeToIndex.at(n3);

        Eigen::Vector3d x = {n1->x(), n2->x(), n3->x()};
        Eigen::Vector3d y = {n1->y(), n2->y(), n3->y()};

        double area = triArea(n1, n2, n3);

        double xc = (x.sum()) / 3.0;
        double yc = (y.sum()) / 3.0;

        double F = f(xc, yc);

        double val = F * area / 3.0;

        b(I1) += val;
        b(I2) += val;
        b(I3) += val;
    }

    return b;
}


// Eigen::VectorXd FemObject::loadVect(const std::list<hed::Edge*>& leading_edges, int np){

//     Eigen::VectorXd b = Eigen::VectorXd::Zero(np);


//     for(const auto& edge : leading_edges){
//         hed::Edge* e1 = edge;
//         hed::Edge* e2 = e1->getNextEdgeInFace();
//         hed::Edge* e3 = e2->getNextEdgeInFace();

//         if (!e1 || !e2 || !e3) continue; // safety check

//         hed::Node* n1 = e1->getSourceNode();
//         hed::Node* n2 = e2->getSourceNode();
//         hed::Node* n3 = e3->getSourceNode();

//         if (!n1 || !n2 || !n3) continue; // skip malformed triangles

//         Eigen::Vector3i loc2glb = {n1->id(),n2->id(),n3->id()};


//         Eigen::Vector3d x_coords = {n1->x(),n2->x(),n3->x()};
//         Eigen::Vector3d y_coords = {n1->y(),n2->y(),n3->y()};
//         double area = triArea(n1,n2,n3);

//         double x_c = (x_coords(0)+x_coords(1)+x_coords(2)) / 3;
//         double y_c = (y_coords(0)+y_coords(1)+y_coords(2)) / 3;

//         double F = f(x_c,y_c);

//         Eigen::Vector3d vec = {1,1,1};
//         Eigen::Vector3d bK = F / 3 * vec * area;

//         for(int j = 0; j < 3; ++j){

//             int J = loc2glb(j);

//             b.coeffRef(J) += bK(j);
//         }


//     }

//     return b;
// }


Eigen::SparseMatrix<double> FemObject::robinMat(
    const std::list<hed::Dart>& boundaries,
    int np,
    const std::unordered_map<const hed::Node*, int>& nodeToIndex)
{
    Eigen::SparseMatrix<double> R(np, np);

    for (auto& boundary : boundaries) {

        hed::Node* n1 = boundary.getNode();
        hed::Node* n2 = boundary.getOppositeNode();

        int I1 = nodeToIndex.at(n1);
        int I2 = nodeToIndex.at(n2);

        double x1 = n1->x(), y1 = n1->y();
        double x2 = n2->x(), y2 = n2->y();

        double length = std::sqrt((x1 - x2)*(x1 - x2) +
                                  (y1 - y2)*(y1 - y2));

        double xc = (x1 + x2) / 2.0;
        double yc = (y1 + y2) / 2.0;

        double k = kappa(xc, yc);

        Eigen::Matrix2d RE;
        RE << 2, 1,
            1, 2;
        RE *= (k * length / 6.0);

        R.coeffRef(I1, I1) += RE(0, 0);
        R.coeffRef(I1, I2) += RE(0, 1);
        R.coeffRef(I2, I1) += RE(1, 0);
        R.coeffRef(I2, I2) += RE(1, 1);
    }

    return R;
}


// Eigen::SparseMatrix<double> FemObject::robinMat(const std::list<hed::Dart>& boundaries, int np){

//     Eigen::SparseMatrix<double> R;
//     R.resize(np,np);

//     for(const auto& boundary: boundaries){
//         hed::Node* n1 = boundary.getNode();
//         hed::Node* n2 = boundary.getOppositeNode();

//         Eigen::Vector2i loc2glb = {n1->id(),n2->id()};

//         Eigen::Vector2d x = {n1->x(),n2->x()};
//         Eigen::Vector2d y = {n1->y(),n2->y()};


//         double length =  std::sqrt((    (  (x(0)  - x(1)) * (x(0)  - x(1)) )    +  (  (y(0)  - y(1)) * (y(0)  - y(1)) ) ));

//         double xc = (x(0) + x(1)) / 2;
//         double yc = (y(0) + y(1)) / 2;


//         double k = kappa(xc,yc);


//         Eigen::Matrix2d mat;
//         mat << 2,1,
//             1,2;

//         Eigen::Matrix2d RE = k / 6 * mat * length;

//         for(int i = 0; i < 2; ++i){
//             for(int j = 0; j < 2; ++j){
//                 // R(loc2glb(i) , loc2glb(j) ) = R(loc2glb(i)  , loc2glb(j)) + RE(i,j);
//                 R.coeffRef(loc2glb(i),loc2glb(j)) += RE(i,j);
//              }
//         }
//     }

//     return R;

// }

Eigen::VectorXd FemObject::robinVect(
    const std::list<hed::Dart>& boundaries,
    int np,
    const std::unordered_map<const hed::Node*, int>& nodeToIndex)
{
    Eigen::VectorXd r = Eigen::VectorXd::Zero(np);

    for (auto& boundary : boundaries) {

        hed::Node* n1 = boundary.getNode();
        hed::Node* n2 = boundary.getOppositeNode();

        int I1 = nodeToIndex.at(n1);
        int I2 = nodeToIndex.at(n2);

        double x1 = n1->x(), y1 = n1->y();
        double x2 = n2->x(), y2 = n2->y();

        double length = std::sqrt((x1 - x2)*(x1 - x2) +
                                  (y1 - y2)*(y1 - y2));

        double xc = (x1 + x2) / 2.0;
        double yc = (y1 + y2) / 2.0;

        double val = (kappa(xc, yc) * gD(xc, yc) + gN(xc, yc)) * length / 2.0;

        r(I1) += val;
        r(I2) += val;
    }

    return r;
}


// Eigen::VectorXd FemObject::robinVect(const std::list<hed::Dart>& boundaries, int np){

//     Eigen::VectorXd r = Eigen::VectorXd::Zero(np);

//     for(const auto& boundary: boundaries){

//         hed::Node* n1 = boundary.getNode();
//         hed::Node* n2 = boundary.getOppositeNode();

//         Eigen::Vector2i loc2glb = {n1->id(),n2->id()};

//         Eigen::Vector2d x = {n1->x(),n2->x()};
//         Eigen::Vector2d y = {n1->y(),n2->y()};


//         double length =  std::sqrt((    (  (x(0)  - x(1)) * (x(0)  - x(1)) )    +  (  (y(0)  - y(1)) * (y(0)  - y(1)) ) ));

//         double xc = (x(0) + x(1)) / 2;
//         double yc = (y(0) + y(1)) / 2;

//         // double tmp = kappa(xc,yc) * gD(xc,yc) - gN(xc,yc);
//         double tmp = kappa(xc,yc) * gD(xc,yc) + gN(xc,yc);

//         Eigen::Vector2d vec = {1.0,1.0};
//         Eigen::Vector2d rE = tmp * vec * length / 2;

//         for(int j = 0; j < 2; j++){



//             r.coeffRef(loc2glb(j)) += rE(j);
//         }

//     }

//     return r;
// }

double FemObject::f(double x, double y){

    ProblemType problemType = this->problemType;

    switch (problemType ) {
    case LAPLACE:
        return 0.0;

    case POISSON:

        return 1.0f;

    case EIGENVALUE:
        return 0.0;

    default:
        return 0.0;
    }

}


double FemObject::gN(double x, double y){
    return 0.0f;
}

double FemObject::gD(double x, double y){

    ProblemType problemType = this->problemType;
    double phi = atan2(y,x);

    switch (problemType) {
    case LAPLACE:
        // phi = atan2(y,x);
        return cos(4*phi);

    case POISSON:
        return y * y * 0.5;

    case HELMHOLTZ:
        return 0.25;

    default:
        return 0;
    }

}

double FemObject::kappa(double x, double y){

    if(this->problemType == HELMHOLTZ || this->problemType == EIGENVALUE){
        if(x > 0.0)
            return 0.0;
        else
            return 1e6;
    }
    return 1e6;
}

void FemObject::setProblemType(ProblemType problemType){
    this->problemType = problemType;
}


void FemObject::solve(){

    const std::list<hed::Node*>* nodes = triang.getNodes();
    if(!nodes || nodes->empty()){

        std::cerr << "[FEMObject::solve] Triangulation has no nodes.\n";

    }

    auto np = static_cast<int>(nodes->size());

    auto& leading_edges = this->triang.getLeadingEdges();

    std::unordered_map<const hed::Node*, int> nodeToIndex;
    nodeToIndex.reserve(np);
    size_t index = 0;

    for (auto node: *nodes){
        nodeToIndex[node] = index++;
    }

    Eigen::SparseMatrix<double> A = stiffMat(leading_edges,np,nodeToIndex);
    Eigen::SparseMatrix<double> M = massMat(leading_edges,np,nodeToIndex);
    Eigen::VectorXd b = loadVect(leading_edges,np, nodeToIndex);

    Eigen::SparseMatrix<double> K;
    Eigen::VectorXd rhs;

    // double lambda = 81; // 81 / 10^2


    double lambda = 81;
    hed::Edge* start = this->triang.getBoundaryEdge();
    if (!start) {
        std::cerr << "No boundary found!\n";
        return;
    }

    hed::Dart bdart(start, true);

    std::list<hed::Dart> boundary;

    ttl::getBoundary(bdart,boundary);

    // std::cout << "Boundary edges = " << boundary.size() << "\n";


    // for (auto& d : boundary) {
    //     auto* nA = d.getNode();
    //     auto* nB = d.getOppositeNode();
    //     std::cout << "Boundary edge: "
    //               << nA->id() << " -> " << nB->id() << "\n";
    // }




    auto R = robinMat(boundary,np, nodeToIndex);

    auto r = robinVect(boundary,np, nodeToIndex);

    // Eigen::SparseMatrix<double> L;




    // for (auto& d : boundary) {
    //     auto* a = d.getNode();
    //     auto* b = d.getOppositeNode();

    //     std::cout << "Boundary edge: (" << a->x() << ", " << a->y() << ") -> ("
    //               << b->x() << ", " << b->y() << ")   kappa = "
    //               << kappa((a->x()+b->x())*0.5, (a->y()+b->y())*0.5) << "\n";
    // }





    switch (problemType) {
    case LAPLACE:

        // FEM formulation:
        // (𝐴 + 𝑅)𝜁 = r
        K = A + R;
        rhs = r;

        break;

    case POISSON:

        // FEM formulation:
        // (𝐴 + 𝑅)𝜁 = b + r

        K = A + R;
        rhs = b + r;

        break;

    case HELMHOLTZ:

        K = A + R - lambda * M;
        rhs = r;

        break;

        // FEM formulation:
        // (A + 𝑅 − 𝜆𝑀)𝜁 = r

    case EIGENVALUE: {
        // Convert to dense
        Eigen::MatrixXd Ld = Eigen::MatrixXd(A + R);
        Eigen::MatrixXd Md = Eigen::MatrixXd(M);

        Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> solver;
        solver.compute(Ld, Md);

        // solver.compute(Md,Ld);

        if (solver.info() != Eigen::Success) {
            std::cerr << "Eigenvalue solve failed\n";
            return;
        }

        Eigen::VectorXd evals = solver.eigenvalues().real();
        Eigen::MatrixXd evecs = solver.eigenvectors().real();

        const int n = static_cast<int>(evals.size());

        // Build a permutation of indices 0..n-1
        std::vector<int> idx(n);
        std::iota(idx.begin(), idx.end(), 0);

        // Sort indices by eigenvalue (ascending)
        std::sort(idx.begin(), idx.end(),
                  [&](int i, int j) { return evals(i) < evals(j); });

        // Optional: print a few eigenvalues to see the spectrum
        // std::cout << "First 10 eigenvalues:\n";
        // for (int k = 0; k < std::min(n, 10); ++k) {
        //     std::cout << k << ": λ = " << evals(idx[k]) << "\n";
        // }

        // Pick the first "physical" mode
        // (skip any zero or tiny/negative eigenvalues if necessary)
        int modeRank = 2;  // 0 = first mode, 1 = second, etc.
        int col = idx[modeRank];

        Eigen::VectorXd zeta = evecs.col(col);

        // Normalize for nicer visualization (optional)
        // double maxAbs = zeta.cwiseAbs().maxCoeff();
        // if (maxAbs > 0.0)
        //     zeta /= maxAbs;

        // Write solution back to nodes using nodeToIndex map
        for (auto* n : *nodes) {
            int id = nodeToIndex.at(n); // use the mapping, NOT n->id()
            double val = zeta(id);
            n->init(n->x(), n->y(), val);
        }
        // double scale = 3.0;  // or 5, 10, etc.
        // for (auto* n : *nodes) {
        //     int id = nodeToIndex.at(n);
        //     double val = scale * zeta(id);
        //     n->init(n->x(), n->y(), val);
        // }

        return;
    }


    // case EIGENVALUE: {
    //     Eigen::MatrixXd Ld = Eigen::MatrixXd(A + R);
    //     Eigen::MatrixXd Md = Eigen::MatrixXd(M);

    //     Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> solver;
    //     solver.compute(Ld,Md);

    //     Eigen::MatrixXcd eigvecC = solver.eigenvectors();
    //     Eigen::MatrixXd eigvec = eigvecC.real();
    //     int modeVal = 1;
    //     if (modeVal >= eigvec.cols()) modeVal = 0;
    //     Eigen::VectorXd zeta = eigvec.col(modeVal);

    //     for (auto* n: *nodes){
    //         int id = nodeToIndex[n];
    //         auto val = zeta(id);
    //         n->init(n->x(),n->y(), val);
    //     }

    //     return;
    // }



        // FEM formulation:
        // (𝐴 + 𝑅) 𝜁 = Λ𝑀𝜁

    }

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;

    solver.compute(K);
    if(solver.info() != Eigen::Success){
        std::cerr << "Factorization failed.\n";
        return;
    }

    Eigen::VectorXd zeta = solver.solve(rhs);

    for(hed::Node* n: *nodes){
        int id = nodeToIndex[n];
        n->init(n->x(),n->y(),zeta(id));
    }


}

double FemObject::u_exact(double x, double y){

    ProblemType problemType = this->problemType;

    switch (problemType) {
    case LAPLACE:{
        auto rho = std::sqrt((x*x) + (y*y)  );
        auto phi = atan2(y,x);
        return std::pow(rho,4.0) * std::cos(4.0*phi);
    }

    case POISSON:{
        return (1 - (x*x)) / 2;
    }

    case HELMHOLTZ:{
        double lambda = 81.0;

        auto cosLambdaX = std::cos(std::sqrt(lambda) * x );
        auto tanLambda = std::tan(std::sqrt(lambda) );
        auto sinLambdaX = std::sin(std::sqrt(lambda) * x );

        return 0.25 * (cosLambdaX + tanLambda * sinLambdaX);
    }

    case EIGENVALUE:
        break;
    }
}

// double FemObject::getError(const std::list<hed::Edge*>& leading_edges,const std::unordered_map<const hed::Node*, int>& nodeToIndex){
double FemObject::getError(){

    // auto& leading_edges = this->triang.getLeadingEdges();


    double error = 0;

    auto& leading_edges = this->triang.getLeadingEdges();

    for (auto edge : leading_edges) {
        hed::Edge* e1 = edge;
        hed::Edge* e2 = e1->getNextEdgeInFace();
        hed::Edge* e3 = e2->getNextEdgeInFace();

        hed::Node* n1 = e1->getSourceNode();
        hed::Node* n2 = e2->getSourceNode();
        hed::Node* n3 = e3->getSourceNode();

        Eigen::Vector3d x = {n1->x(), n2->x(), n3->x()};
        Eigen::Vector3d y = {n1->y(), n2->y(), n3->y()};

        double area = triArea(n1, n2, n3);

        double xc = (x.sum()) / 3.0;
        double yc = (y.sum()) / 3.0;

        auto ubar = u_exact(xc,yc);
        auto uhbar =  (  n1->z() + n2->z() + n3->z() )    / 3.0 ;

        auto eK = ( ubar - uhbar   ) * ( ubar - uhbar   ) * area;

        error += eK;
    }
    return std::sqrt(error);

}

int FemObject::getDoF(){
    return this->triang.getNodes()->size();
}


std::string FemObject::problemTypeToString(ProblemType pt) {
    switch (pt) {
    case LAPLACE:    return "laplace";
    case POISSON:    return "poisson";
    case HELMHOLTZ:  return "helmholtz";
    case EIGENVALUE: return "eigenvalue";
    default:         return "unknown";
    }
}


