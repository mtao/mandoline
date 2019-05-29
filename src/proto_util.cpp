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
}
