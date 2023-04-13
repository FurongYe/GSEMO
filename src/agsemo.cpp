#include "agsemo.hpp"

/***
 * ./main problem_id dimension lambda p budget runs algorithm_name
 */
int main(int argc, char *argv[])
{
    if (argc < 6)
    {
        std::cerr << "Some parameters are missing" << std::endl;
    }

    int problem_id = std::atoi(argv[1]);
    int dimension = std::atoi(argv[2]);
    int lambda = std::atoi(argv[3]);
    double pm = std::atoi(argv[4]) / static_cast<double>(dimension);
    int budget = std::atoi(argv[5]);
    int runs = std::atoi(argv[6]);
    std::string algorithm_name = argv[7];

    // int problem_id = 4;
    // int dimension = 100;
    // int lambda = 1;
    // double pm = 1 / static_cast<double>(dimension);
    // int budget = 10000000;
    // int runs = 10;
    // std::string algorithm_name = "AGSEMO";
   
    using namespace ioh::common;
    random::seed(10);

    bool force_flip = true;

    auto problem = get_problem(problem_id, dimension);
    
    
    
    std::string exp_name = algorithm_name + "L" + argv[4] + "P" + argv[5];
    
    auto logger =
        ioh::logger::Analyzer({ioh::trigger::always}, {ioh::watch::violation}, "MO", exp_name, exp_name, exp_name);

    problem->attach_logger(logger);

    for (int i = 0; i < runs; ++i)
    {
        AGSEMO<OptimizationType::MAX, OptimizationType::MAX> opt(budget, force_flip, pm, lambda, algorithm_name);
        auto p = opt(problem);
        print(p);
        std::cout << p.size() << std::endl;
        std::cout << *problem << std::endl;
        problem->reset();
    }
}
