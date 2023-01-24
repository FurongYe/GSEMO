#include "gsemo.hpp"

/***
 * ./main problem_id dimension static/TwoRate/varctrl/logNormal lambda p adapt_metric budget runs
 * adapt_metric  1: hypervolume 2: IGD 3: NOPareto
 */
int main(int argc, char *argv[])
{
    if (argc < 6)
    {
        std::cerr << "Some parameters are missing" << std::endl;
    }

    int problem_id = std::atoi(argv[1]);
    int dimension = std::atoi(argv[2]);
    std::string algorithm_name = argv[3];
    int lambda = std::atoi(argv[4]);
    double pm = std::atoi(argv[5]) / static_cast<double>(dimension);
    int adapt_metric = std::atoi(argv[6]);
    int budget = std::atoi(argv[7]);
    int runs = std::atoi(argv[8]);

    // int problem_id = 1;
    // int dimension = 100;
    // int lambda = 2;
    // double pm = 1 / static_cast<double>(dimension);
    // int budget = 100000;
    // int runs = 2;

    using namespace ioh::common;
    random::seed(10);

    bool force_flip = true;

    auto problem = get_problem(problem_id, dimension);
    
    // std::string algorithm_name = "TwoRate";
    std::string metric = (adapt_metric == 1 ? "HV" : (adapt_metric == 2 ? "IGD" : "NUM"));
    std::string exp_name = algorithm_name + "L" + argv[4] + "P" + argv[5] + metric;
    // std::string exp_name = "test";
    auto logger =
        ioh::logger::Analyzer({ioh::trigger::always}, {ioh::watch::violation}, "MO", exp_name, exp_name, exp_name);

    problem->attach_logger(logger);

    for (int i = 0; i < runs; ++i)
    {
        GSEMO<OptimizationType::MAX, OptimizationType::MAX> opt(budget, force_flip, pm, lambda, adapt_metric, algorithm_name);
        auto p = opt(problem);
        print(p);
        std::cout << p.size() << std::endl;
        std::cout << *problem << std::endl;
        problem->reset();
    }
}
