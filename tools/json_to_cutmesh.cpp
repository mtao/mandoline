#include <nlohmann/json.hpp>
#include "utils/json_to_cutmesh.hpp"






int main(int argc, char * argv[]) {




    auto ccm = json_to_cutmesh(std::string(argv[1]));
    ccm.write(argv[2]);

    return 0;

}
