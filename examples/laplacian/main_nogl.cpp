#include "mandoline/mesh3.hpp"
#include "setup.h"

int main(int argc, char* argv[]) {
    auto ccm = read(argv[1]);
    std::cout << divergence(ccm).transpose() << std::endl;
    return 0;


}
