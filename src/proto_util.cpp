#include "mandoline/proto_util.hpp"
namespace mandoline::protobuf {
    void serialize(const mtao::Vec2i& a, CutMeshProto::Vec2i& b) {
        b.set_x(a.x());
        b.set_y(a.y());
    }
    void serialize(const mtao::Vec3d& a, CutMeshProto::Vec3d& b) {
        b.set_x(a.x());
        b.set_y(a.y());
        b.set_z(a.z());
    }
    void serialize(const mtao::Vec3i& a, CutMeshProto::Vec3i& b) {
        b.set_x(a.x());
        b.set_y(a.y());
        b.set_z(a.z());
    }
    void serialize(const std::array<int,2>& a, CutMeshProto::Vec2i& b) {
        b.set_x(a[0]);
        b.set_y(a[1]);
    }
    void serialize(const std::array<int,3>& a, CutMeshProto::Vec3i& b) {
        b.set_x(a[0]);
        b.set_y(a[1]);
        b.set_z(a[2]);
    }
    void serialize(const std::array<double,3>& a, CutMeshProto::Vec3d& b) {
        b.set_x(a[0]);
        b.set_y(a[1]);
        b.set_z(a[2]);
    }
    void deserialize(const CutMeshProto::Vec2i& a, std::array<int,2>& b   ) {
        b[0] = a.x();
        b[1] = a.y();
    }
    void deserialize(const CutMeshProto::Vec3i& a, std::array<int,3>& b   ) {
        b[0] = a.x();
        b[1] = a.y();
        b[2] = a.z();
    }
    void deserialize(const CutMeshProto::Vec3d& a, std::array<double,3>& b) {
        b[0] = a.x();
        b[1] = a.y();
        b[2] = a.z();
    }
    mtao::Vec3d deserialize(const CutMeshProto::Vec3d& a) {
        mtao::Vec3d b;
        b.x() = a.x();
        b.y() = a.y();
        b.z() = a.z();
        return b;
    }
    mtao::Vec3i deserialize(const CutMeshProto::Vec3i& a) {
        mtao::Vec3i b;
        b.x() = a.x();
        b.y() = a.y();
        b.z() = a.z();
        return b;
    }
    mtao::Vec2i deserialize(const CutMeshProto::Vec2i& a) {
        mtao::Vec2i b;
        b.x() = a.x();
        b.y() = a.y();
        return b;
    }

    void deserialize(const CutMeshProto::Vertex& a, Vertex<3>& b   ) {
        b.coord[0] = a.i();
        b.coord[1] = a.j();
        b.coord[2] = a.k();
        b.quot(0) = a.u();
        b.quot(1) = a.v();
        b.quot(2) = a.w();
        b.clamped_indices[0] = a.ci();
        b.clamped_indices[1] = a.cj();
        b.clamped_indices[2] = a.ck();
    }
    void serialize(const Vertex<3>& a, CutMeshProto::Vertex& b) {
        b.set_i(a.coord[0]);
        b.set_j(a.coord[1]);
        b.set_k(a.coord[2]);
        b.set_u(a.quot(0));
        b.set_v(a.quot(1));
        b.set_w(a.quot(2));
        b.set_ci(a.clamped_indices[0]);
        b.set_cj(a.clamped_indices[1]);
        b.set_ck(a.clamped_indices[2]);
    }
}
