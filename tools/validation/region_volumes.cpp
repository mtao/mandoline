#define __gmp_const const
#include "cutmesh_validation.hpp"
#include <igl/copyleft/cgal/extract_cells.h>
using namespace mandoline;
std::array<int,2> pcwn_count(const CutCellMesh<3>& ccm) {

    return {{is_pcwn,is_not_pcwn}};
}
bool pcwn_check(const CutCellMesh<3>& ccm) {
    return std::get<1>(pcwn_count(ccm)) == 0;
}
bool volume_check(const mandoline::CutCellMesh<3>&) {
}
