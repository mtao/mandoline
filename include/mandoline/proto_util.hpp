#pragma once
#include <mtao/types.hpp>
#include "mandoline/vertex.hpp"
#include "cutmesh.pb.h"

namespace mandoline::protobuf {
    void serialize(const mtao::Vec3i& a, CutMeshProto::Vec3i& b);
    void serialize(const mtao::Vec3d& a, CutMeshProto::Vec3d& b);
    void serialize(const std::array<int,3>& a, CutMeshProto::Vec3i& b);
    void serialize(const std::array<double,3>& a, CutMeshProto::Vec3d& b);
    void serialize(const mtao::Vec2i& a, CutMeshProto::Vec2i& b);
    void serialize(const std::array<int,2>& a, CutMeshProto::Vec2i& b);
    void serialize(const Vertex<3>& a, CutMeshProto::Vertex& b);

    mtao::Vec3d deserialize(const CutMeshProto::Vec3d& a);
    mtao::Vec3i deserialize(const CutMeshProto::Vec3i& a);
    mtao::Vec2i deserialize(const CutMeshProto::Vec2i& a);

    void deserialize(const CutMeshProto::Vec3i& a, std::array<int,3>& b   );
    void deserialize(const CutMeshProto::Vec3d& a, std::array<double,3>& b);
    void deserialize(const CutMeshProto::Vec2i& a, std::array<int,2>& b   );
    void deserialize(const CutMeshProto::Vertex& a, Vertex<3>& b   );
}
