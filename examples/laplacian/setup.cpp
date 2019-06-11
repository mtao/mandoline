#include "setup.h"


mandoline::CutCellMesh<3> read(const std::string& filename) {
    return mandoline::CutCellMesh<3>::from_proto(filename);
}


std::map<int,int> face_regions(const mandoline::CutCellMesh<3>& ccm) {
    std::map<int,int> FR;
    for(int i = 0; i < ccm.faces().size(); ++i) {
        FR[i] = 0;
    }
    for(auto&& c: ccm.cells()){
        const int region = c.region;
        for(auto&& [f,s]: c) {
            int& r2 = FR[f];
            r2 = std::max(r2,region);
        }
    }
    return FR;
}


mtao::VecXd flux(const mandoline::CutCellMesh<3>& ccm) {
    auto FR = face_regions(ccm);
    mtao::VecXd FV = ccm.face_volumes();
    auto r = ccm.regions();
    int rcount = *std::max_element(r.begin(),r.end());
    std::vector<double> RV(rcount+1,0);
    std::cout << "RCount: " << rcount << std::endl;
    std::copy(RV.begin(),RV.end(),std::ostream_iterator<int>(std::cout,", "));
    std::cout << std::endl;
    for(auto&& [i,f]: mtao::iterator::enumerate(ccm.faces())) {
        //std::cout << FR.at(i) << "+=" << FV(i) << std::endl;;
        //std::cout << RV[FR.at(i)] << " => ";
        RV[FR.at(i)] += FV(i);
        //std::cout << RV[FR.at(i)] << std::endl;;
        //if(f.is_axial_face()) {
        //    FV(i) = 0;
        //}
    }
    std::copy(RV.begin(),RV.end(),std::ostream_iterator<int>(std::cout,", "));
    std::cout << std::endl;
    for(auto&& [i,f]: mtao::iterator::enumerate(ccm.faces())) {
        if(RV[FR.at(i)] > 0) {
        FV(i) /= RV[FR.at(i)];
        }
    }
    //std::cout << FV.transpose() << std::endl;
    std::cout << "Flux rows: " << FV.rows() << std::endl;
    return FV;
}
mtao::VecXd divergence(const mandoline::CutCellMesh<3>& ccm) {
    return boundary(ccm).transpose() * flux(ccm);
}
Eigen::SparseMatrix<double> boundary(const mandoline::CutCellMesh<3>& ccm) {
    auto B = ccm.boundary();
    std::cout << "Boundary size: " << B.rows() << "," << B.cols() << std::endl;
    return B;
}
Eigen::SparseMatrix<double> laplacian(const mandoline::CutCellMesh<3>& ccm) {
    return {};
}
