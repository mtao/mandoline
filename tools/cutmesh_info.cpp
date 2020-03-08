#include <mtao/cmdline_parser.hpp>
#include <mandoline/tools/cutmesh_info.hpp>
#include <iostream>


int main(int argc, char * argv[]) {
    mtao::logging::make_logger().set_level(mtao::logging::Level::Off);
    if(argc < 1)
    {
        std::cout << "cutmesh_info <filename1> <filename2> ..." << std::endl;
    }
    for(int i = 1; i < argc; ++i)
    {
        if(i > 1)
        {
            std::cout << std::endl << std::endl;
        }
        mandoline::tools::print_file_info(argv[i]);
    }


    return 0;

}
