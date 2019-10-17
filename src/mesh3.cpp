#include "mandoline/mesh3.hpp"
#include <mtao/iterator/enumerate.hpp>
#include <mtao/geometry/bounding_box.hpp>
#include <mtao/geometry/volume.h>
#include <mtao/logging/logger.hpp>
#include <mtao/geometry/mesh/halfedge.hpp>
#include <mtao/geometry/prune_vertices.hpp>
#include <mtao/geometry/mesh/separate_elements.hpp>
#include <mtao/geometry/kdtree.hpp>
#include <mtao/iterator/range.hpp>
#include <igl/winding_number.h>
#include "mandoline/diffgeo_utils.hpp"
#include "mandoline/proto_util.hpp"


namespace mandoline {

    size_t CutCellMesh<3>::cell_size() const {
        return m_adaptive_grid.num_cells() + m_cells.size();
    }
    size_t CutCellMesh<3>::cut_face_size() const {
        return m_faces.size();
    }
    size_t CutCellMesh<3>::face_size() const {
        return m_adaptive_grid.num_faces() + cut_face_size();
    }
    size_t CutCellMesh<3>::cut_cell_size() const { 
        return m_cells.size(); 
    }
    bool   CutCellMesh<3>::is_cut_cell(int index) const {
        return index >= 0 && index < m_cells.size(); 
    }
    bool   CutCellMesh<3>::is_adaptive_cell(int index) const { 
        return index >= 0 && index >= m_cells.size(); 
    }
    bool CutCellMesh<3>::is_mesh_face(int idx) const {
        return m_mesh_faces.find(idx) != m_mesh_faces.end();
    }
    auto CutCellMesh<3>::cell_volumes() const -> VecX{
        /*
           VecX V(cell_size());
           V.topRows(StaggeredGrid::cell_size()).array() = dx().prod();
           */
        VecX V(cells().size());
        V.setZero();
        auto Vs = vertices();
        mtao::VecXd face_brep_vols(m_faces.size());
        for(auto&& [i,f]: mtao::iterator::enumerate(m_faces)) {
            face_brep_vols(i) = f.brep_volume(Vs);
        }

        for(auto&& [k,c]: mtao::iterator::enumerate(m_cells)) {
            V(k) = c.volume(face_brep_vols);
            //V(k) = c.volume(Vs,faces);
            //V(k) = mtao::geometry::brep_volume(Vs,c.triangulated(faces));
            //std::cout << (V(k) / std::abs(c.volume(Vs,faces))) << std::endl;
        }


        return mtao::eigen::vstack(V,adaptive_grid().cell_volumes());
    }
    auto CutCellMesh<3>::face_centroids() const -> mtao::ColVecs3d{
        /*
           VecX V(cell_size());
           V.topRows(StaggeredGrid::cell_size()).array() = dx().prod();
           */
        mtao::ColVecs3d ret(3,m_faces.size());
        auto Vs = vertices();
        for(auto&& [i,f]: mtao::iterator::enumerate(m_faces)) {
            auto v = ret.col(i);
            v.setZero();
            int count = 0;
            for(auto&& ind: f.indices) {
                for(auto&& j: ind) {
                    v += Vs.col(j);
                    count++;
                }
            }
            if(count > 0) {
                v /= count;
            }
        }
        auto r = adaptive_grid().boundary_centroids();
        return mtao::eigen::hstack(ret,r);
    }
    auto CutCellMesh<3>::cell_centroids() const -> mtao::ColVecs3d{
        /*
           VecX V(cell_size());
           V.topRows(StaggeredGrid::cell_size()).array() = dx().prod();
           */
        mtao::ColVecs3d face_brep_cents(3,face_size());
        auto Vs = vertices();
        for(auto&& [i,f]: mtao::iterator::enumerate(m_faces)) {
            face_brep_cents.col(i) = f.brep_centroid(Vs);
            //face_brep_cents.col(i) = f.brep_volume(Vs) * f.brep_centroid(Vs);
        }
        mtao::ColVecs3d V(3,cell_size());
        V.setZero();
        auto vols = cell_volumes();
        for(auto&& [k,c]: mtao::iterator::enumerate(m_cells)) {
            V.col(k) = c.moment(face_brep_cents);
            //V.col(k) = c.moment(face_brep_cents) / vols(k);
        }


        adaptive_grid().cell_centroids(V);
        return V;
    }
    auto CutCellMesh<3>::dual_vertices() const -> ColVecs {
        return cell_centroids();
    }
    int CutCellMesh<3>::local_grid_cell_index(const VecCRef& p) const {
        auto [c,q] = StaggeredGrid::coord(p);
        if(!StaggeredGrid::cell_grid().valid_index(c)) {
            return -1;
        }

        return -1;
    }
    int CutCellMesh<3>::world_grid_cell_index(const VecCRef& p) const {
        return local_grid_cell_index(vertex_grid().local_coord(p));
    }
    std::vector<std::set<int>> CutCellMesh<3>::cell_faces() const {
        std::vector<std::set<int>> ret(m_cells.size());
        std::transform(m_cells.begin(),m_cells.end(), ret.begin(), [](auto&& V) {
                std::set<int> ret;
                std::transform(V.begin(),V.end(), std::inserter(ret,ret.end()), [](auto&& v) {
                        return std::get<0>(v);
                        });

                return ret;
                });
        return ret;

    }
    std::set<int> CutCellMesh<3>::cell_faces(int index) const {
        auto& V = m_cells.at(index);
        std::set<int> ret;
        std::transform(V.begin(),V.end(), std::inserter(ret,ret.end()), [](auto&& v) {
                return std::get<0>(v);
                });

        return ret;

    }

    void CutCellMesh<3>::write(const std::string& prefix) const {
        GOOGLE_PROTOBUF_VERIFY_VERSION;
        std::stringstream ss;
        auto CS = StaggeredGrid::cell_shape();
        {
            ss << prefix << "." << CS[0] << "_" << CS[1] << "_" << CS[2];
            ss << ".cutmesh";
        }
        //std::cout << "Output filename: " << ss.str() << std::endl;
        std::ofstream ofs(ss.str(), std::ios::binary);
        CutMeshProto cmp;
        serialize(cmp);
        cmp.SerializeToOstream(&ofs);


    }

    void CutCellMesh<3>::write_obj_regions(const std::string& filename) const {
        std::vector<int> R = regions();
        std::map<int,std::set<int>> Rs;
        for(auto&& [i,r]: mtao::iterator::enumerate(R)) {
            Rs[r].insert(i);
        }
        for(auto&& [i,inds]: Rs) {
            write_obj(filename + "_regions",inds,i,false);
        }

    }
    void CutCellMesh<3>::write_obj_flaps(const std::string& filename) const {
        std::set<int> inds;
        for(int i = 0; i < cell_size(); ++i) {
            inds.insert(i);
        }
        write_obj(filename,inds,{},false,false,true);
    }
    void CutCellMesh<3>::write_obj(const std::string& filename) const {
        std::set<int> inds;
        for(int i = 0; i < cell_size(); ++i) {
            inds.insert(i);
        }
        write_obj(filename,inds,{},false);
    }
    void CutCellMesh<3>::write_mesh_obj_separate(const std::string& prefix) const {
        auto R = regions();
        int i = 0;
#pragma omp parallel for
        for(i = 0; i < cell_size(); ++i) {
            write_obj(prefix,{i},R[i], true,false,false,true);
        }
    }
    void CutCellMesh<3>::write_obj_separate(const std::string& prefix, bool separate_flaps) const {
        auto R = regions();
        int i = 0;
#pragma omp parallel for
        for(i = 0; i < cell_size(); ++i) {
            //write_obj(prefix,i);
            if(separate_flaps) {
                write_obj(prefix,{i},R[i], true, true,false);
                write_obj(prefix,{i},R[i], true, false,true);
            } else {
                std::cout << "Writing " << i << std::endl;
                write_obj(prefix,{i},R[i], true, true,true);
            }
        }
    }
    void CutCellMesh<3>::write_obj(const std::string& prefix, const std::set<int>& indices, const std::optional<int>& region, bool show_indices, bool show_base, bool show_flaps, bool m_mesh_faces) const {
        std::stringstream ss;
        auto CS = StaggeredGrid::cell_shape();

        {
            ss << prefix;
            if(m_mesh_faces) {
                ss << "-mesh";
                if(!(show_base && show_flaps)) {
                    if(show_base) {
                        ss << "-base";
                    } else if(show_flaps) {
                        ss << "-flap";
                    }
                }
            }
            if(show_indices) {
                ss <<  "_" ;
                for(auto&& index: indices) {
                    ss << index << "_" ;
                    if(is_cut_cell(index)) {
                        auto c = m_cells.at(index).grid_cell;
                        ss << "(" << c[0] << "_" << c[1] << "_" << c[2] << ")";
                    }
                }
            } else {
                ss <<  "_" ;
            }
            if(region) {
                ss << "r" << *region<< ".obj";
            } else {
                ss << ".obj";
            }
        }

        auto V = vertices();

        std::vector<mtao::ColVecs3i> Fs(indices.size());
        std::vector<mtao::ColVecs3d> Vs(indices.size());
        int vertex_offset = V.cols();
        for(auto&& [V,F,index]: mtao::iterator::zip(Vs,Fs,indices)) {
            if(m_mesh_faces) {
                if(is_cut_cell(index)) {
                    std::vector<mtao::ColVecs3i> FF;
                    for(auto&& [fidx,s]: m_cells.at(index)) {
                        auto&& f = m_faces[fidx];
                        if(f.is_mesh_face()) {
                            FF.push_back(*f.triangulation);
                        }
                    }
                    F = mtao::eigen::hstack_iter(FF.begin(),FF.end());
                }
            } else {
                if(is_cut_cell(index)) {

                    std::tie(V,F) = triangulated_cell(index,show_base,show_flaps);
                } else {
                    if(show_base) {
                        F = m_adaptive_grid.triangulated(index);
                    }
                }
            }
            if(is_cut_cell(index)) {
                auto c = m_cells.at(index).grid_cell;
                std::cout << "Cell: " << c[0] << " " << c[1] << " " << c[2] << std::endl;
            }
        }


        //std::cout << "Output filename: " << ss.str() << std::endl;
        {
            bool empty = true;
            for(auto&& f: Fs) {
                if(f.size() > 0) {
                    empty = false;
                }
            }
            if(empty) {
                std::cout << "Empty cell!" << std::endl;
                return;
            }
        }

        std::ofstream ofs(ss.str());
        std::cout << "Output file: " << ss.str() << std::endl;
        //if(normalize_output) {
        //    auto bb = mtao::geometry::bounding_box(V);

        //    V.colwise() -= origin();
        //    V.array().colwise() /= dx().array().cast<double>();
        //} else if(normalize_unit_grid_cell) {
        //}

        if(region) {
            ofs << "# region " << *region << std::endl;
        }



        std::tie(V,Fs) = mtao::geometry::mesh::compactify(V,Fs);
        for(int i = 0; i < V.cols(); ++i) {
            auto v = V.col(i);
            ofs << "v " << v.transpose() << std::endl;
        }

        for(auto&& [F,index]: mtao::iterator::zip(Fs,indices)) {
            if(is_cut_cell(index)) {
                ofs << "# cell" << index << ": ";
                for(auto&& [c,b]: m_cells.at(index)) {
                    ofs << c << " ";
                }
            } else {
                ofs << "# adaptive cube" << index << ": ";
                auto&& c = m_adaptive_grid.cell(index);
                ofs <<  mtao::eigen::stl2eigen(c.corner()).transpose() << "  :: " << c.width();
            }
            ofs << std::endl;

            for(int i = 0; i < F.cols(); ++i) {
                auto f = F.col(i).array()+1;
                ofs << "f " << f.transpose() << std::endl;
            }
        }
    }

    std::vector<bool> CutCellMesh<3>::boundary_faces() const {
        std::vector<bool> ret(m_faces.size(),false);
        for(auto&& [m,_]: m_mesh_faces) {
            ret[m] = true;
        }
        return ret;
    }

    std::vector<int> CutCellMesh<3>::regions(bool boundary_sign_regions) const {
        std::vector<int> ret(cell_size(),1);

        if(boundary_sign_regions) {
            std::cout << "Using boundary sign regions!" << std::endl;
            for(auto&& [idx,c]: mtao::iterator::enumerate(m_cells)) {
                Edge counts{{0,0}};// 1 -1 
                for(auto&& [f,b]: c) {
                    if(is_mesh_face(f)) {
                        counts[b?0:1]++;
                    }
                }
                if(counts[0] + counts[1] > 0) {
                    if(counts[0] > counts[1]) {
                        ret[idx] = 2;
                    } else {
                        ret[idx] = 3;
                    }
                }
            }
        } else {
            std::cout << "Baking regions!" << std::endl;
            for(int i = 0; i < m_cells.size(); ++i) {
                ret[i] = m_cells[i].region;
            }
            for(auto&& [c,r]: m_adaptive_grid_regions) {
                ret[c] = r;
            }
        }
        return ret;

    }

    std::vector<std::array<std::set<int>,2>> CutCellMesh<3>::face_regions() const {
        auto R = regions();
        std::vector<std::array<std::set<int>,2>> ret(R.size());
        for(auto&& c: cells()) {
            for(auto&& [fidx,s]: c) {
                auto& f = faces()[fidx];
                if(f.is_mesh_face()) {
                    ret[R[c.index]][s?0:1].insert(fidx);
                }
            }
        }
        return ret;
    }
    std::vector<std::array<std::set<int>,2>> CutCellMesh<3>::orig_face_regions() const {
        auto R = face_regions();
        std::vector<std::array<std::set<int>,2>> ret(R.size());
        for(auto&& [Fsp,OFsp]: mtao::iterator::zip(R,ret)) {
            for(auto&& [Fs,OFs]: mtao::iterator::zip(Fsp,OFsp)) {
                for(auto&& f: Fs) {
                    assert(faces()[f].is_mesh_face());
                    OFs.insert(faces()[f].as_face_id());
                }
            }
        }
        return ret;
    }

    mtao::ColVecs3d CutCellMesh<3>::region_centroids() const {
        //delay dividing by 3 until the very end...
        auto R = orig_face_regions();
        mtao::ColVecs3d Cs(3,origF().cols());
        mtao::VecXd vols(origF().cols());
        for(int i = 0; i < Cs.cols(); ++i) {
            auto f = origF().col(i);
            auto c = Cs.col(i);
            mtao::Mat3d A;
            for(int j = 0; j < 3; ++j) {
                A.col(j) = origV().col(f(j));
            }
            vols(i) = A.determinant();
            c = vols(i)  * A.rowwise().sum();
        }

        mtao::ColVecs3d C(3,R.size());
        mtao::VecXd V(R.size());
        C.setZero();
        V.setZero();
        for(auto&& [rid, rp]: mtao::iterator::enumerate(R)) {
            auto c = C.col(rid);
            auto&& v = V(rid);
            for(auto&& nf: rp[1]) {
                c -= Cs.col(nf);
                v -= vols(nf);
            }
            for(auto&& pf: rp[0]) {
                c += Cs.col(pf);
                v += vols(pf);
            }
        }


        //TODO: handle outside boundaries of the grid using the adaptive grid!
        auto bb = bbox();
        double gv = bb.sizes().prod();
        mtao::Vec3d gc = 1.5 * (bb.min()+bb.min()) * gv;
        C.col(0) += gc;
        V(0) += gv;

        for(int i = 0; i < C.cols(); ++i) {
            C.col(i) /= 3 * V(i);
        }
        return C;
    }


    void CutCellMesh<3>::write_mesh_surface_obj(const std::string& prefix) const {
        std::stringstream ss;
        auto CS = StaggeredGrid::cell_shape();

        {
            ss << prefix;
            ss << "-mesh";
            ss << prefix << "." << CS[0] << "_" << CS[1] << "_" << CS[2];
            ss << ".obj";
        }
        std::vector<mtao::ColVecs3i> Fs;
        for(auto&& f: m_faces) {
            if(f.is_mesh_face()) {
                if(f.triangulation) {
                    Fs.push_back(*f.triangulation);
                }
            }
        }
        std::ofstream ofs(ss.str());
        auto V = vertices();
        std::tie(V,Fs) = mtao::geometry::mesh::compactify(V,Fs);
        for(int i = 0; i < V.cols(); ++i) {
            auto v = V.col(i);
            ofs << "v " << v.transpose() << std::endl;
        }

        for(auto&& F: Fs) {
            for(int i = 0; i < F.cols(); ++i) {
                auto f = F.col(i).array()+1;
                ofs << "f " << f.transpose() << std::endl;
            }
        }
    }

    CutCellMesh<3> CutCellMesh<3>::from_file(const std::string& filename) {
        return from_proto(filename);
    }


    void  CutCellMesh<3>::serialize(CutMeshProto& cmp) const {
        for(int i = 0; i < cut_vertex_size(); ++i) {
            protobuf::serialize(cut_vertex(i),*cmp.add_vertices());
        }
        for(int i = 0; i < cut_edge_size(); ++i) {
            protobuf::serialize(m_cut_edges.col(i),*cmp.add_edges());
        }
        for(int i = 0; i < m_origV.cols(); ++i) {
            protobuf::serialize(m_origV.col(i),*cmp.add_origv());
        }
        for(int i = 0; i < m_origF.cols(); ++i) {
            protobuf::serialize(m_origF.col(i),*cmp.add_origf());
        }
        protobuf::serialize(origin(),*cmp.mutable_origin());
        protobuf::serialize(dx(),*cmp.mutable_dx());
        protobuf::serialize(shape(),*cmp.mutable_shape());
        for(auto&& f: m_faces) {
            f.serialize(*cmp.add_faces());
        }
        for(auto&& c: m_cells) {
            c.serialize(*cmp.add_cells());
        }
        for(auto&& i: m_folded_faces) {
            cmp.add_foldedfaces(i);
        }
        auto&& mf = *cmp.mutable_mesh_faces();
        for(auto&& [idx, bmf]: m_mesh_faces) {
            auto& b = mf[idx];
            b.set_parent_id(bmf.parent_fid);
            for(int j = 0; j < bmf.barys.cols(); ++j) {
                protobuf::serialize(bmf.barys.col(j),*b.add_barycentric_coordinates());

            }

        }
        {
            auto& cmap = *cmp.mutable_cubes();
            for(auto&& [c,cell]: m_adaptive_grid.cells()) {
                cell.serialize(cmap[c]);
            }
            auto& rmap = *cmp.mutable_cube_regions();
            for(auto&& [a,b]: m_adaptive_grid_regions) {
                rmap[a] = b;
            }
        }
    }
    CutCellMesh<3> CutCellMesh<3>::from_proto(const std::string& filename) {
        std::ifstream ifs(filename,std::ios::binary);
        CutMeshProto cmp;
        cmp.ParseFromIstream(&ifs);
        return from_proto(cmp);
    }
    CutCellMesh<3> CutCellMesh<3>::from_proto(const CutMeshProto& cmp) {

        mtao::Vec3d o,dx;
        std::array<int,3> s;
        o = protobuf::deserialize(cmp.origin());
        dx = protobuf::deserialize(cmp.dx());
        protobuf::deserialize(cmp.shape(),s);
        CutCellMesh<3> ret = CutCellMesh<3>::StaggeredGrid(s,dx,o);
        //ret.m_face_volumes.resize(20);

        ret.m_cut_vertices.resize(cmp.vertices().size());
        for(int i = 0; i < ret.cut_vertex_size(); ++i) {
            protobuf::deserialize(cmp.vertices(i),ret.m_cut_vertices[i]);
        }
        ret.m_cut_edges.resize(2,cmp.edges().size());
        for(int i = 0; i < ret.cut_edge_size(); ++i) {
            ret.m_cut_edges.col(i) = protobuf::deserialize(cmp.edges(i));
        }
        ret.m_origV.resize(3,cmp.origv().size());
        for(int i = 0; i < ret.m_origV.cols(); ++i) {
            ret.m_origV.col(i) = protobuf::deserialize(cmp.origv(i));
        }
        ret.m_origF.resize(3,cmp.origf().size());
        for(int i = 0; i < ret.m_origF.cols(); ++i) {
            ret.m_origF.col(i) = protobuf::deserialize(cmp.origf(i));
        }
        ret.m_faces.resize(cmp.faces().size());
        for(int i = 0; i < cmp.faces().size(); ++i) {
            ret.m_faces[i] = CutFace<3>::from_proto(cmp.faces(i));
            ret.m_faces[i].update_mask(ret.cut_vertices(),ret.vertex_grid());
        }
        ret.m_cells.resize(cmp.cells().size());
        for(int i = 0; i < cmp.cells().size(); ++i) {
            ret.m_cells[i] = CutCell::from_proto(cmp.cells(i));
        }
        std::copy(cmp.foldedfaces().begin(),cmp.foldedfaces().end(),std::inserter(ret.m_folded_faces,ret.m_folded_faces.end()));


        for(auto&& [idx,btf]: cmp.mesh_faces()) {
            int pid = btf.parent_id();
            size_t bssize = btf.barycentric_coordinates_size();
            mtao::ColVecs3d B(3,bssize);
            for(int i = 0; i < bssize; ++i) {
                B.col(i) = protobuf::deserialize(btf.barycentric_coordinates(i));
            }
            ret.m_mesh_faces[idx] = {B,pid};
        }

        for(auto&& [a,b]: cmp.cubes()) {
            ret.m_adaptive_grid.m_cells[a] = AdaptiveGrid::Cell::from_proto(b);
        }
        ret.m_adaptive_grid.make_faces();
        for(auto&& [a,b]: cmp.cube_regions()) {
            ret.m_adaptive_grid_regions[a] = b;
        }

        return ret;
    }
    std::array<mtao::ColVecs2d,3> CutCellMesh<3>::compute_subVs() const {
        auto V = vertices();
        std::array<mtao::ColVecs2d,3> subVs;
        for(int d = 0; d < 3; ++d) {
            int n0 = (d+1)%3;
            int n1 = (d+2)%3;

            auto&& subV = subVs[d];
            subV.resize(2,V.cols());
            subV.row(0) = V.row(n0);
            subV.row(1) = V.row(n1);

        }
        return subVs;
    }
    std::tuple<mtao::ColVecs3d,mtao::ColVecs3i> CutCellMesh<3>::triangulate_face(int face_index) const {
        mtao::logging::warn() << "Inefficient use of triangulation!  try caching your triangulations";
        std::array<mtao::ColVecs2d,3> subVs = compute_subVs();
        auto [V,F] = m_faces[face_index].triangulate(subVs,false);
        if(V.cols() > 0) {
            return {vertex_grid().world_coord(V),F};
        } else {
            return {{},F};
        }
    }

    void CutCellMesh<3>::triangulate_faces() {
        std::array<mtao::ColVecs2d,3> subVs = compute_subVs();

        int i = 0;
#pragma omp parallel for
        for(i = 0; i < m_faces.size(); ++i) {
            auto&& face = m_faces[i];
            face.cache_triangulation(subVs,true);
        }
    }

    mtao::VecXd CutCellMesh<3>::dual_edge_lengths() const {
        VecX DL = VecX::Zero(face_size());

        auto& dx = Base::dx();
        auto g = adaptive_grid().grid();
        for(auto&& c: cells()) {
            auto& gc = c.grid_cell;
            for(auto&& [fidx,s]: c) {
                auto& f = faces()[fidx];
                if(f.external_boundary) {
                    auto [cid,s] = *f.external_boundary;
                    auto&& c = adaptive_grid().cell(g.get(cid));
                    mtao::Vec3d C = c.center();
                    C.array() -= .5;
                    C -= mtao::eigen::stl2eigen(gc).cast<double>();

                    DL(fidx) = (dx.asDiagonal() * C).norm();
                }
            }
        }
        for(auto&& [fidx,f]: mtao::iterator::enumerate(faces())) {
            if(f.is_axial_face() && !f.external_boundary) {
                int ba = f.as_axial_axis();
                DL(fidx) = dx(ba);
            }
        }
        auto ADL = adaptive_grid().dual_edge_lengths();
        DL.tail(ADL.rows()) = ADL;
        return DL;
    }

    auto CutCellMesh<3>::face_volumes(bool from_triangulation) const -> VecX{
        VecX FV(faces().size());
        if(from_triangulation) {
            for(auto&& [i,face]: mtao::iterator::enumerate(faces())) {
                auto V = vertices();
                if(face.triangulation) {
                    auto&& T = *face.triangulation;
                    FV(i) = mtao::geometry::volumes(V,T).sum();
                }
            }
            FV = mtao::eigen::vstack(FV,adaptive_grid().face_volumes());
        } else {
            //Use barycentric for tri-mesh cut-faces, planar areas for axial cut-faces, and let the adaptive grid do it's thing
            auto trimesh_vols = mtao::geometry::volumes(m_origV,m_origF);

            if(trimesh_vols.size() > 0) {
                FV = face_barycentric_volume_matrix() * trimesh_vols;
                auto subVs = compute_subVs();
                for(auto&& [i,f]: mtao::iterator::enumerate(faces())) {
                    if(f.is_axial_face()) {
                        auto [dim,coord] = f.as_axial_id();
                        auto& vol = FV(i) = 0;
                        auto& V = subVs[dim];
                        for(auto&& indices: f.indices) {
                            vol += mtao::geometry::curve_volume(V,indices);
                        }
                        if(vol < 0) {
                            mtao::logging::warn() << "Negative volume! warn mtao because this shouldn't happen";
                            vol = -vol;
                        }

                    }
                    if(!std::isfinite(FV(i))) {
                        FV(i) = 0;
                    }
                }
            } else {
                FV.resize(adaptive_grid().num_faces());
            }
            FV.tail(adaptive_grid().num_faces()) = adaptive_grid().face_volumes();
        }

        return FV;
    }

    mtao::VecXd CutCellMesh<3>::mesh_face_mask() const {
        mtao::VecXd M = mtao::VecXd::Ones(face_size());
        for(auto&& [i,f]: mtao::iterator::enumerate(faces())) {
            if(f.is_mesh_face()) {
                M(i) = 0;
            }
        }
        return M;
    }
    mtao::VecXd CutCellMesh<3>::primal_hodge2() const {
        auto PV = face_volumes();
        auto DV = dual_edge_lengths();
        mtao::VecXd CV = (PV.array() > 1e-5).select(DV.cwiseQuotient(PV),0);
        for(int i = 0; i < CV.size(); ++i) {
            if(!std::isfinite(CV(i))) {
                CV(i) = 0;
            }
        }
        return CV;
    }

    mtao::VecXd CutCellMesh<3>::dual_hodge2() const {
        auto PV = face_volumes();
        auto DV = dual_edge_lengths();
        mtao::VecXd  CV = (DV.array() > 1e-5).select(PV.cwiseQuotient(DV),0);
        for(int i = 0; i < CV.size(); ++i) {
            if(!std::isfinite(CV(i))) {
                CV(i) = 0;
            }
        }
        return CV;
    }
    mtao::VecXd CutCellMesh<3>::dual_hodge3() const {
        auto CV = cell_volumes();
        for(int i = 0; i < CV.size(); ++i) {
            CV(i) = (std::abs(CV(i)) < 1e-5) ? 0 : (1. / CV(i));
            if(!std::isfinite(CV(i))) {
                CV(i) = 0;
            }
        }
        return CV;
    }
    mtao::VecXd CutCellMesh<3>::primal_hodge3() const {
        auto CV = cell_volumes();
        for(int i = 0; i < CV.size(); ++i) {
            if(!std::isfinite(CV(i))) {
                CV(i) = 0;
            }
        }
        return CV;
    }

    Eigen::SparseMatrix<double> CutCellMesh<3>::trilinear_matrix() const {
        Eigen::SparseMatrix<double> A(vertex_size(), StaggeredGrid::vertex_size());
        std::vector<Eigen::Triplet<double>> trips;
        trips.reserve(StaggeredGrid::vertex_size() + 8 * cut_vertex_size());
        for(int i = 0; i < StaggeredGrid::vertex_size(); ++i) {
            trips.emplace_back(i,i,1);
        }
        int offset = StaggeredGrid::vertex_size();
        for(auto&& [idx,v]: mtao::iterator::enumerate(cut_vertices())) {
            int index = idx + offset;
            coord_type a = v.coord;
            for(int i = 0; i < 2; ++i) {
                a[0] = v.coord[0] + i;
                double vi = i==0?(1-v.quot(0)):v.quot(0);
                for(int j = 0; j < 2; ++j) {
                    a[1] = v.coord[1] + j;
                    double vij = vi * (j==0?(1-v.quot(1)):v.quot(1));
                    for(int k = 0; k < 2; ++k) {
                        a[2] = v.coord[2] + k;
                        double vijk = vij * (k==0?(1-v.quot(2)):v.quot(2));
                        trips.emplace_back(index,vertex_index(a),vijk);
                    }
                }
            }
        }
        A.setFromTriplets(trips.begin(),trips.end());
        return A;
    }
    Eigen::SparseMatrix<double> CutCellMesh<3>::face_grid_volume_matrix() const {
        auto trips = m_adaptive_grid.grid_face_projection(m_faces.size());
        auto FV = face_volumes();
        mtao::Vec3d gfv;
        gfv(0) = dx()(1) * dx()(2);
        gfv(1) = dx()(0) * dx()(2);
        gfv(2) = dx()(0) * dx()(1);
        Eigen::SparseMatrix<double> A(this->face_size(),this->form_size<2>());
        for(auto&& t: trips) {
            const int row = t.row();
            const int col = t.col();
        }
        for(auto&& [i,face]: mtao::iterator::enumerate(faces())) {
            if(face.count() == 1) {
                int axis = face.bound_axis();
                constexpr static int maxval = std::numeric_limits<int>::max();
                coord_type c{{maxval,maxval,maxval}};
                for(auto&& ind: face.indices) {
                    for(auto&& i: ind) {
                        auto v = masked_vertex(i).coord;
                        for(auto&& [a,b]: mtao::iterator::zip(c,v)) {
                            a = std::min(a,b);
                        }
                    }
                }
                const int row = i;
                const int col = staggered_index<2>(c,axis);
                double value = FV(i) / gfv(axis) * face.N(axis);//face.N(axis) should be a unit vector either facing up or down....
                trips.emplace_back(row,col,value);
            }
        }
        A.setFromTriplets(trips.begin(),trips.end());
        //mtao::VecXd sums = A * mtao::VecXd::Zero(A.cols());
        //sums = (sums.array().abs() > 1e-10).select(1.0 / sums.array(), 0);
        //A = sums.asDiagonal() * A;
        return A;
    }

    Eigen::SparseMatrix<double> CutCellMesh<3>::barycentric_matrix() const {
        Eigen::SparseMatrix<double> A(vertex_size(), m_origV.cols());
        std::vector<Eigen::Triplet<double>> trips;
        std::map<std::array<int,2>,double> mp;
        for(auto&& [fid, btf]: m_mesh_faces) {
            auto t = btf.sparse_entries(m_faces[fid], m_origF);
            std::copy(t.begin(),t.end(),std::inserter(mp,mp.end()));
        }
        trips.reserve(mp.size());
        for(auto&& [pr,v]: mp) {
            auto [a,b] = pr;
            trips.emplace_back(a,b,v);
        }
        A.setFromTriplets(trips.begin(),trips.end());
        mtao::VecXd sums = A * mtao::VecXd::Ones(A.cols());
        sums = (sums.array().abs() > 1e-10).select(1.0 / sums.array(), 0);
        A = sums.asDiagonal() * A;
        return A;
    }
    Eigen::SparseMatrix<double> CutCellMesh<3>::face_barycentric_volume_matrix() const {
        int face_size = 0;
        //artifact from before i passed in m_origF
        if(m_origF.size() == 0) {
            for(auto&& f: faces()) {
                if(f.is_mesh_face()) {
                    face_size = std::max<int>(face_size,f.as_face_id());
                }
            }
            face_size++;
        } else {
            face_size = m_origF.cols();
        }
        Eigen::SparseMatrix<double> A(this->face_size(),face_size);
        std::vector<Eigen::Triplet<double>> trips;
        for(auto&& [fid, btf]: m_mesh_faces) {
            double vol = btf.volume() * 2;//proportion of 2 is required because barycentric coordinates live in a unit triangle
            trips.emplace_back(fid, btf.parent_fid, vol);
        }
        A.setFromTriplets(trips.begin(),trips.end());
        //for(int i = 0; i < A.cols(); ++i) {
        //    A.col(i) /= A.col(i).sum();
        //}
        return A;
    }
    Eigen::SparseMatrix<double> CutCellMesh<3>::boundary() const {
        auto trips = m_adaptive_grid.boundary_triplets(m_faces.size());
        Eigen::SparseMatrix<double> B(face_size(),cell_size());

        auto g = adaptive_grid().grid();

        for(auto&& c: cells()) {
            int region = c.region;
            for(auto&& [fidx,s]: c) {
                auto& f = faces()[fidx];
                if(f.is_axial_face()) {
                    int a = cell_shape()[f.as_axial_axis()];
                    int v = f.as_axial_coord();
                    if(v > 0 && v < a) {
                        trips.emplace_back(fidx,c.index, s?1:-1);
                    }
                } else {
                    trips.emplace_back(fidx,c.index, s?1:-1);
                }
            }
        }
        for(auto&& [fidx,f]: mtao::iterator::enumerate(faces())) {
            if(f.external_boundary) {
                auto [c,s] = *f.external_boundary;
                trips.emplace_back(fidx,g.get(c), s?-1:1);
            }
        }


        B.setFromTriplets(trips.begin(),trips.end());
        return B;
    }
    auto CutCellMesh<3>::active_cell_mask() const -> GridDatab {
        auto mask = GridDatab::Constant(true,StaggeredGrid::cell_shape());
        auto C = cell_centroids();
        for(auto&& [i,cell]: mtao::iterator::enumerate(m_cells)) {
            auto cent = C.col(i);
            auto [c,q] = coord(cent);
            mask(c) = false;
        }
        return mask;
    }
    auto CutCellMesh<3>::active_cells() const -> std::set<coord_type> {
        auto C = cell_centroids();
        std::set<coord_type> ret;
        for(auto&& [i,cell]: mtao::iterator::enumerate(m_cells)) {
            auto cent = C.col(i);
            auto [c,q] = coord(cent);
            ret.insert(c);
        }
        return ret;
    }
    size_t CutCellMesh<3>::active_cell_count() const {
        return active_cells().size();
    }

    std::tuple<mtao::ColVecs3d, mtao::ColVecs3i> CutCellMesh<3>::triangulated_cell(int idx, bool use_base, bool use_flap) const {
        if(is_cut_cell(idx)) {
            std::vector<mtao::ColVecs3d> mVs;
            std::vector<mtao::ColVecs3i> mFs;
            for(auto&& [fidx,b]: m_cells[idx]) {
                bool is_flap = m_folded_faces.find(fidx) != m_folded_faces.end();
                bool is_base = !is_flap;
                if((use_flap && is_flap) || (use_base && !is_flap)) {
                    auto& tri = m_faces[fidx].triangulation;
                    mtao::ColVecs3d V;
                    mtao::ColVecs3i F;
                    if(tri) {
                        F = *tri;
                        auto& vert = m_faces[fidx].triangulated_vertices;
                        if(vert) {
                            V = *vert;
                        }
                    } else {
                        std::tie(V,F) = triangulate_face(fidx);
                    }
                    if(!b) {
                        auto tmp = F.row(0).eval();
                        F.row(0) = F.row(1);
                        F.row(1) = tmp;
                    }
                    mVs.emplace_back(std::move(V));
                    mFs.emplace_back(std::move(F));
                }
            }
            return {mtao::eigen::hstack_iter(mVs.begin(),mVs.end()),mtao::eigen::hstack_iter(mFs.begin(),mFs.end())};
            /*
            if(mFs.size() > 0) {
                    return {mtao::eigen::hstack_iter(mVs.begin(),mVs.end()),mtao::eigen::hstack_iter(mFs.begin(),mFs.end())};
                if(mVs.size() > 0) {
                    return {mtao::eigen::hstack_iter(mVs.begin(),mVs.end()),mtao::eigen::hstack_iter(mFs.begin(),mFs.end())};
                } else {
                    return {{},mtao::eigen::hstack_iter(mFs.begin(),mFs.end())};
                }
            }
            */
        } else {
            if(use_base) {
                return {{},m_adaptive_grid.triangulated(idx)};
            }
        }
        return {};
    }

    std::tuple<mtao::ColVecs3d,mtao::ColVecs3i> CutCellMesh<3>::compact_triangulated_cell(int cell_index) const {
        auto [V,F] = triangulated_cell(cell_index,true,true);
        return mtao::geometry::mesh::compactify(mtao::eigen::hstack(vertices(),V),F);
    }
    std::tuple<mtao::ColVecs3d,mtao::ColVecs3i> CutCellMesh<3>::compact_triangulated_face(int face_index) const {
        auto& f = m_faces[face_index];
        if(!f.triangulation) {
            mtao::logging::warn() << "Triangulate mesh first!";
            return {};
        } else {
            if(f.triangulated_vertices) {
                return mtao::geometry::mesh::compactify(mtao::eigen::hstack(vertices(),vertex_grid().world_coord(*f.triangulated_vertices)),*f.triangulation);
            } else {
                return mtao::geometry::mesh::compactify(vertices(),*f.triangulation);
            }
        }
    }

    auto CutCellMesh<3>::cells_by_grid_cell() const -> std::map<coord_type,std::set<int>> {
        std::map<coord_type,std::set<int>> R;
        for(auto&& [i,c]: mtao::iterator::enumerate(cells())) {
            R[c.grid_cell].insert(i);
        }
        return R;
    }
    std::set<int> CutCellMesh<3>::cells_in_grid_cell(const coord_type& c) const {
        std::set<int> R;
        for(auto&& [i,cell]: mtao::iterator::enumerate(cells())) {
            if(cell.grid_cell == c) {
            R.insert(i);
            }
        }
        return R;
    }
    int CutCellMesh<3>::get_cell_index(const VecCRef& p) const {
        auto v = local_coord(p);
        //check if its  in an adaptive grid cell, tehn we can just use that cell
        if(int ret = m_adaptive_grid.get_cell_index(v); ret != -1) {
            return ret;
        } else {

            auto [c,q] = vertex_grid().coord(p);
            auto cell_indices = cells_in_grid_cell(c);
            for(auto&& ci: cell_indices) {
                auto&& cell = cells()(ci);
                if(cell.is_inside(p)) {
                    return ci;
                }
            }
            //TODO
            //if I haven't returned yet then either this ccm is bad
            //or i'm too close to an edge. lets not assume that for now
        }
        return -1;
    }
}
