#pragma once
#include <balsa/eigen/types.hpp>
#include "mandoline/vertex.hpp"
#include "cutmesh.pb.h"

namespace mandoline::protobuf {
void serialize(const balsa::eigen::Vec3i &a, Vec3i &b);
void serialize(const balsa::eigen::Vec3d &a, Vec3d &b);
void serialize(const std::array<int, 3> &a, Vec3i &b);
void serialize(const std::array<double, 3> &a, Vec3d &b);
void serialize(const balsa::eigen::Vec2i &a, Vec2i &b);
void serialize(const std::array<int, 2> &a, Vec2i &b);
void serialize(const mandoline::Vertex<3> &a, Vertex &b);

balsa::eigen::Vec3d deserialize(const Vec3d &a);
balsa::eigen::Vec3i deserialize(const Vec3i &a);
balsa::eigen::Vec2i deserialize(const Vec2i &a);

void deserialize(const Vec3i &a, std::array<int, 3> &b);
void deserialize(const Vec3d &a, std::array<double, 3> &b);
void deserialize(const Vec2i &a, std::array<int, 2> &b);
void deserialize(const Vertex &a, mandoline::Vertex<3> &b);
}// namespace mandoline::protobuf
