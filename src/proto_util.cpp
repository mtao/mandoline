#include "mandoline/proto_util.hpp"
namespace mandoline::protobuf {
void serialize(const balsa::eigen::Vec2i &a, Vec2i &b) {
    b.set_x(a.x());
    b.set_y(a.y());
}
void serialize(const balsa::eigen::Vec3d &a, Vec3d &b) {
    b.set_x(a.x());
    b.set_y(a.y());
    b.set_z(a.z());
}
void serialize(const balsa::eigen::Vec3i &a, Vec3i &b) {
    b.set_x(a.x());
    b.set_y(a.y());
    b.set_z(a.z());
}
void serialize(const std::array<int, 2> &a, Vec2i &b) {
    b.set_x(a[0]);
    b.set_y(a[1]);
}
void serialize(const std::array<int, 3> &a, Vec3i &b) {
    b.set_x(a[0]);
    b.set_y(a[1]);
    b.set_z(a[2]);
}
void serialize(const std::array<double, 3> &a, Vec3d &b) {
    b.set_x(a[0]);
    b.set_y(a[1]);
    b.set_z(a[2]);
}
void deserialize(const Vec2i &a, std::array<int, 2> &b) {
    b[0] = a.x();
    b[1] = a.y();
}
void deserialize(const Vec3i &a, std::array<int, 3> &b) {
    b[0] = a.x();
    b[1] = a.y();
    b[2] = a.z();
}
void deserialize(const Vec3d &a, std::array<double, 3> &b) {
    b[0] = a.x();
    b[1] = a.y();
    b[2] = a.z();
}
balsa::eigen::Vec3d deserialize(const Vec3d &a) {
    balsa::eigen::Vec3d b;
    b.x() = a.x();
    b.y() = a.y();
    b.z() = a.z();
    return b;
}
balsa::eigen::Vec3i deserialize(const Vec3i &a) {
    balsa::eigen::Vec3i b;
    b.x() = a.x();
    b.y() = a.y();
    b.z() = a.z();
    return b;
}
balsa::eigen::Vec2i deserialize(const Vec2i &a) {
    balsa::eigen::Vec2i b;
    b.x() = a.x();
    b.y() = a.y();
    return b;
}

void deserialize(const Vertex &a, mandoline::Vertex<3> &b) {
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
void serialize(const mandoline::Vertex<3> &a, Vertex &b) {
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
}// namespace mandoline::protobuf
