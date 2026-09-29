#include "patch_query_validation.hpp"
#include <cstdlib>
#include <exception>
#include <iostream>
int main()
{
    try
    {
        rainbow::tests::test_bilinear_queries();
        rainbow::tests::test_query_collector();
        rainbow::tests::test_direction_grid();
        std::cout<<"Patch query CPU: analytic, 12000 independent-reference queries, collector and angular grid passed.\n";
        return EXIT_SUCCESS;
    }
    catch(const std::exception& e){std::cerr<<"Error: "<<e.what()<<'\n';return EXIT_FAILURE;}
}
