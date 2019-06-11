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
        return m_adaptive_grid.cells.size() + m_cells.size();
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
        VecX V(cell_size());
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
        double scale = dx().array().prod();
        for(auto&& [i,c]:  m_adaptive_grid.cells) {
            V(i) =  scale * std::pow<double>(c.width(),3);
        }


        return V;
    }
    auto CutCellMesh<3>::cell_centroids() const -> mtao::ColVecs3d{
        /*
           VecX V(cell_size());
           V.topRows(StaggeredGrid::cell_size()).array() = dx().prod();
           */
        mtao::ColVecs3d face_brep_cents(3,m_faces.size());
        auto Vs = vertices();
        for(auto&& [i,f]: mtao::iterator::enumerate(m_faces)) {
            face_brep_cents.col(i) = f.brep_centroid(Vs);
            //face_brep_cents.col(i) = f.brep_volume(Vs) * f.brep_centroid(Vs);
        }
        mtao::ColVecs3d V(3,m_cells.size());
        V.setZero();
        auto vols = cell_volumes();
        for(auto&& [k,c]: mtao::iterator::enumerate(m_cells)) {
            V.col(k) = c.moment(face_brep_cents);
            //V.col(k) = c.moment(face_brep_cents) / vols(k);
        }


        return V;
    }
    auto CutCellMesh<3>::dual_vertices() const -> ColVecs {
        return cell_centroids();
    }
    int CutCellMesh<3>::cell_index(const VecCRef& p) const {
        auto [c,q] = StaggeredGrid::coord(p);
        if(!StaggeredGrid::cell_grid().valid_index(c)) {
            return -1;
        }

        return -1;
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

        std::vector<mtao::ColVecs3i> Fs(indices.size());
        for(auto&& [F,index]: mtao::iterator::zip(Fs,indices)) {
            if(m_mesh_faces) {
                if(is_cut_cell(index)) {
                    std::cout << "cell: " << index << std::endl;
                    std::vector<mtao::ColVecs3i> FF;
                    for(auto&& [fidx,s]: m_cells.at(index)) {
                        std::cout << fidx << " ";
                        auto&& f = m_faces[fidx];
                        if(f.is_mesh_face()) {
                            FF.push_back(*f.triangulation);
                        }
                    }
                    std::cout << std::endl;
                    F = mtao::eigen::hstack_iter(FF.begin(),FF.end());
                }
            } else {
                if(is_cut_cell(index)) {
                    F = triangulated_cell(index,show_base,show_flaps);
                } else {
                    if(show_base) {
                        F = m_adaptive_grid.triangulated(index);
                    }
                }
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
                //std::cout << "Empty cell!" << std::endl;
                return;
            }
        }

        std::ofstream ofs(ss.str());
        auto V = vertices();
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
                auto&& c = m_adaptive_grid.cells.at(index);
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
                    std::cout << f << "/" << m_mesh_faces.size() << std::endl;
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
        for(int i = 0; i < cut_vertices.cols(); ++i) {
            protobuf::serialize(cut_vertices.col(i),*cmp.add_vertices());
        }
        for(int i = 0; i < cut_edges.cols(); ++i) {
            protobuf::serialize(cut_edges.col(i),*cmp.add_edges());
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
            for(auto&& [c,cell]: m_adaptive_grid.cells) {
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

        ret.cut_vertices.resize(3,cmp.vertices().size());
        for(int i = 0; i < ret.cut_vertices.cols(); ++i) {
            ret.cut_vertices.col(i) = protobuf::deserialize(cmp.vertices(i));
        }
        ret.cut_edges.resize(2,cmp.edges().size());
        for(int i = 0; i < ret.cut_edges.cols(); ++i) {
            ret.cut_edges.col(i) = protobuf::deserialize(cmp.edges(i));
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
            ret.m_adaptive_grid.cells[a] = AdaptiveGrid::Cell::from_proto(b);
        }
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
    mtao::ColVecs3i CutCellMesh<3>::triangulate_face(int face_index) const {
        mtao::logging::warn() << "Inefficient use of triangulation!  try caching your triangulations";
        std::array<mtao::ColVecs2d,3> subVs = compute_subVs();
        return m_faces[face_index].triangulate(subVs);
    }

    void CutCellMesh<3>::triangulate_faces() {
        std::array<mtao::ColVecs2d,3> subVs = compute_subVs();

        int i = 0;
#pragma omp parallel for
        for(i = 0; i < m_faces.size(); ++i) {
            auto&& face = m_faces[i];
            face.cache_triangulation(subVs);
        }
    }

    auto CutCellMesh<3>::face_volumes(bool from_triangulation) const -> VecX{
        VecX ret(faces().size());
        if(from_triangulation) {
            for(auto&& [i,face]: mtao::iterator::enumerate(faces())) {
                auto V = vertices();
                if(face.triangulation) {
                    auto&& T = *face.triangulation;
                    ret(i) = mtao::geometry::volumes(V,T).sum();
                }
            }
        } else {
            auto trimesh_vols = mtao::geometry::volumes(m_origV,m_origF);

            ret = face_barycentric_volume_matrix() * trimesh_vols;
            auto subVs = compute_subVs();
            for(auto&& [i,f]: mtao::iterator::enumerate(faces())) {
                if(!std::isfinite(ret(i))) {
                    ret(i) = 0;
                }
                if(f.is_axial_face()) {
                    auto [dim,coord] = f.as_axial_id();
                    auto& vol = ret(i) = 0;
                    auto& V = subVs[dim];
                    for(auto&& indices: f.indices) {
                        vol += mtao::geometry::curve_volume(V,indices);
                    }
                    if(vol < 0) {
                        mtao::logging::warn() << "Negative volume! warn mtao because this shouldn't happen";
                        vol = -vol;
                    }

                }
            }
        }
        return ret;
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
        return A;
    }
    Eigen::SparseMatrix<double> CutCellMesh<3>::face_barycentric_volume_matrix() const {
        int face_size = 0;
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
        Eigen::SparseMatrix<double> A(m_faces.size(),face_size);
        std::vector<Eigen::Triplet<double>> trips;
        for(auto&& [fid, btf]: m_mesh_faces) {
            double vol = btf.volume() * 2;//proportion of 2 is required because barycentric coordinates live in a unit triangle
            trips.emplace_back(fid, btf.parent_fid, vol);
        }
        A.setFromTriplets(trips.begin(),trips.end());
        return A;
    }
    Eigen::SparseMatrix<double> CutCellMesh<3>::boundary() const {
        auto trips = m_adaptive_grid.boundary_triplets(m_faces.size());
        int num_faces = faces().size();
        int num_cells = cells().size();
        for(auto&& t: trips) {
            num_faces = std::max<int>(t.row()+1,num_faces);
            num_cells = std::max<int>(t.col()+1,num_faces);
        }
        Eigen::SparseMatrix<double> B(num_faces,num_cells);


        for(auto&& c: cells()) {
            int region = c.region;
            for(auto&& [fidx,s]: c) {
                trips.emplace_back(fidx,c.index, s?1:-1);
            }
        }
        for(auto&& [fidx,f]: mtao::iterator::enumerate(faces())) {
            if(f.external_boundary) {
                auto [c,s] = *f.external_boundary;
                trips.emplace_back(fidx,c, s?1:-1);
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

    mtao::ColVecs3i CutCellMesh<3>::triangulated_cell(int idx, bool use_base, bool use_flap) const {
        if(is_cut_cell(idx)) {
            std::vector<mtao::ColVecs3i> mFs;
            for(auto&& [fidx,b]: m_cells[idx]) {
                bool is_flap = m_folded_faces.find(fidx) != m_folded_faces.end();
                bool is_base = !is_flap;
                if((use_flap && is_flap) || (use_base && !is_flap)) {
                    auto& tri = m_faces[fidx].triangulation;
                    mtao::ColVecs3i F = tri?(*tri).eval():triangulate_face(fidx);
                    if(!b) {
                        auto tmp = F.row(0).eval();
                        F.row(0) = F.row(1);
                        F.row(1) = tmp;
                    }
                    mFs.emplace_back(std::move(F));
                }
            }
            if(mFs.size() > 0) {
                return mtao::eigen::hstack_iter(mFs.begin(),mFs.end());
            }
        } else {
            if(use_base) {
                return m_adaptive_grid.triangulated(idx);
            }
        }
        return {};
    }

    std::tuple<mtao::ColVecs3d,mtao::ColVecs3i> CutCellMesh<3>::compact_triangulated_cell(int cell_index) const {
        auto F = triangulated_cell(cell_index,true,true);
        return mtao::geometry::mesh::compactify(vertices(),F);
    }
    std::tuple<mtao::ColVecs3d,mtao::ColVecs3i> CutCellMesh<3>::compact_triangulated_face(int face_index) const {
        auto& f = m_faces[face_index];
        if(!f.triangulation) {
            mtao::logging::warn() << "Triangulate mesh first!";
            return {};
        } else {
            return mtao::geometry::mesh::compactify(vertices(),*f.triangulation);
        }
    }
}
