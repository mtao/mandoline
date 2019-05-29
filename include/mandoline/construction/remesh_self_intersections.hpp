#include <mtao/types.hpp>

namespace mandoline::construction {
    std::tuple<mtao::ColVecs3d,mtao::ColVecs3i> remesh_self_intersections(const mtao::ColVecs3d& V,const mtao::ColVecs3i& F);
}
