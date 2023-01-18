#include "ioh.hpp"
#include "problems_gecco23.hpp"

int dimension;

namespace operators
{
    template <ioh::common::OptimizationType>
    bool is_worse(const double a, const double b);

    template <>
    bool is_worse<ioh::common::OptimizationType::MIN>(const double a, const double b)
    {
        return a >= b;
    }

    template <>
    bool is_worse<ioh::common::OptimizationType::MAX>(const double a, const double b)
    {
        return a <= b;
    }

    template <ioh::common::OptimizationType>
    bool is_better(const double a, const double b);

    template <>
    bool is_better<ioh::common::OptimizationType::MIN>(const double a, const double b)
    {
        return a <= b;
    }

    template <>
    bool is_better<ioh::common::OptimizationType::MAX>(const double a, const double b)
    {
        return a >= b;
    }

    template <ioh::common::OptimizationType>
    double distance_(const double a, const double b);

    template <>
    double distance_<ioh::common::OptimizationType::MAX>(const double a, const double b)
    {
        return a - b;
    }

    template <>
    double distance_<ioh::common::OptimizationType::MIN>(const double a, const double b)
    {
        return a - b;
    }


    /// \fn sampleNFromM
    /// \brief Sampling n different indexes from length maninclude
    void sampleNFromM(std::vector<size_t> &sampled_number, size_t n, size_t m)
    {
        if (sampled_number.size() != 0)
        {
            sampled_number.clear();
        }

        if (n == 0)
        {
            std::clog << "sampled zero number" << std::endl;
        }

        size_t randPos;
        sampled_number.reserve(n);

        if (n > m / 2)
        { /// If n is larger than m/2, we sample random indexes by reordering a permutation.
            std::vector<size_t> population;
            population.reserve(m);
            for (size_t i = 0; i < m; ++i)
            {
                population.push_back(i);
            }

            int temp;
            for (size_t i = m - 1; i > 0; --i)
            {
                randPos = ioh::common::random::integer(0, i);
                temp = population[i];
                population[i] = population[randPos];
                population[randPos] = temp;
                sampled_number.push_back(population[i]);
                if (m - i - 1 == n - 1)
                {
                    break;
                }
            }
            if (n == m)
            {
                sampled_number.push_back(population[0]);
            }
        }
        else
        { /// If n is smaller than m/2, we sample indexes repeatly until getting different values.
            bool resample = false;
            for (size_t i = 0; i != n; ++i)
            {
                do
                {
                    resample = false;
                    randPos = ioh::common::random::integer(0, m - 1);
                    for (size_t j = 0; j != i; ++j)
                    {
                        if (randPos == sampled_number[j])
                        {
                            resample = true;
                            break;
                        }
                    }
                } while (resample);
                sampled_number.push_back(randPos);
            }
        }
    };
} // namespace operators

template <ioh::common::OptimizationType... Args>
struct MultiSolution
{
    static constexpr int S = sizeof...(Args);
    std::array<double, S> y;
    std::vector<int> x;

    bool operator==(const MultiSolution<Args...> &other) const
    {
        for (size_t i = 0; i < S; i++)
            if (y[i] != other.y[i])
                return false;
        return true;
    }

    bool dominated_by(const MultiSolution<Args...> &other) const
    {
        bool dominated = true;
        int i = 0;
        (
            [&] {
                dominated &= operators::is_worse<Args>(y[i], other.y[i]);
                i++;
            }(),
            ...);
        return dominated;
    }

    bool strictly_dominates(const MultiSolution<Args...> &other) const
    {
        bool dominates = true;
        int i = 0;
        (
            [&] {
                dominates &= operators::is_better<Args>(y[i], other.y[i]);
                i++;
            }(),
            ...);
        return dominates and !(*this == other);
    }

    double distance(const MultiSolution<Args...> &other) const
    {
        double dis;
        int i = 0;
        (
            [&] {
                dis += operators::distance_<Args>(y[i], other.y[i]);
                i++;
            }(),
            ...);
        return dis;
    }

    void print() const
    {
        std::cout << "(";
        for (const auto &fi : y)
            std::cout << fi << ",";
        std::cout << ")";
    }
};

template <ioh::common::OptimizationType... Args>
void print(const std::vector<MultiSolution<Args...>> &p)
{
    std::cout << "{";
    for (const auto &z : p)
        z.print();
    std::cout << "}\n";
}

void bitflip(std::vector<int> &x, const double pm)
{
    for (size_t i = 0; i < x.size(); i++)
    {
        const auto ri = ioh::common::random::real();
        if (ri < pm)
        {
            x[i] = abs(x[i] - 1);
        }
    }
}

void bitflip_index(std::vector<int> &x, const int flip_bits_number)
{
    std::vector<size_t> flip_index;
    operators::sampleNFromM(flip_index, flip_bits_number, x.size());
    for (const auto &pos : flip_index)
        x[pos] = abs(x[pos] - 1);
}

void bitflip_binom(std::vector<int> &x, const double pm, const int m = 0)
{
    std::binomial_distribution<> d(x.size(), pm);
    const auto n = std::max(m, d(ioh::common::random::gen));
    for (const auto &pos : ioh::common::random::integers(n, 0, (int)x.size() - 1))
        x[pos] = abs(x[pos] - 1);
}


template <ioh::common::OptimizationType... Args>
struct GSEMO
{
    using GSolution = MultiSolution<Args...>;

    int budget;
    bool force_flip;
    double pm_star;
    int lambda;
    int verbose_rate;
    std::string algorithm_name;

    double r = 1;
    double pm;
    double F = 0.98;
    int c = 0;

    std::normal_distribution<> standard_norm{0, 1};
    // Standard mersenne_twister_engine seeded with rd()

    GSEMO(const int b = 6000, const bool ff = true, const double pm = 0.01, const int lambda = 1,
          const std::string algorithm_name = "static", const int v = 0) :
        budget(b),
        pm_star(pm), lambda(lambda), algorithm_name(algorithm_name), force_flip(ff), verbose_rate(v)
    {
    }


    std::vector<GSolution> operator()(const std::shared_ptr<ioh::problem::Integer> &problem)
    {
        auto x = std::vector<int>(problem->meta_data().n_variables, 0.0);
        r = pm_star * problem->meta_data().n_variables;
        bitflip(x, 0.5);
        std::vector<GSolution> p{eval(x, problem)};
        --budget;
        while (budget > 0)
        {
            std::vector<GSolution> pnew;
            std::vector<double> pmNew;
            for (int i = 1; i <= lambda; i++)
            {
                auto ri = ioh::common::random::integer(0, p.size() - 1);
                auto xi = p[ri].x;

                /***
                 * Pm relization
                 */
                int l;
                if (algorithm_name == "static")
                {
                    pm = pm_star;
                    std::binomial_distribution<> d{problem->meta_data().n_variables, pm};
                    l = d(ioh::common::random::gen);
                    while (l < 1 ){
                        l = d(ioh::common::random::gen);
                    }
                }
                else if (algorithm_name == "tworate")
                {
                    if (i < lambda / 2)
                    {
                        pm = r / 2.0 / problem->meta_data().n_variables;
                    }
                    else
                    {
                        pm = r * 2.0 / problem->meta_data().n_variables;
                    }
                    std::binomial_distribution<> d{problem->meta_data().n_variables, pm};
                    
                    l = d(ioh::common::random::gen);
                    while (l < 1){
                        l = d(ioh::common::random::gen);
                    }
                    pmNew.push_back(pm);
                }

                else if (algorithm_name == "varctrl")
                {
                    std::normal_distribution<> d{
                        r, r * (1 - r / static_cast<double>(problem->meta_data().n_variables)) * pow(F, c)};
                    l = std::round(d(ioh::common::random::gen));
                    while (l < 1 || l > problem->meta_data().n_variables / 2)
                    {
                        l = std::round(d(ioh::common::random::gen));
                    }
                    pmNew.push_back(l);
                }
                else if (algorithm_name == "logNormal")
                {
                    pm = 1.0 /
                        (1.0 + (((1.0 - pm_star) / pm_star) * exp(0.22 * standard_norm(ioh::common::random::gen))));
                    std::binomial_distribution<> d{problem->meta_data().n_variables, pm};
                    l = d(ioh::common::random::gen);
                    while (l < 1) {
                    l = d(ioh::common::random::gen);
                    }
                    pmNew.push_back(pm);
                }


                bitflip_index(xi, l);
                auto si = eval(xi, problem);
                if (--budget <= 0)
                    break;

                if (problem->state().optimum_found)
                    break;


                // std::vector<GSolution> pnew{si};

                pnew.push_back(si);

                if (verbose_rate and i % verbose_rate == 0)
                    print(p);
            }

            /***
             * For self-adaptation
             */
            
            if (algorithm_name != "static")
            {
                std::map<double, int> success_time;
                std::map<double, double> distances;
                auto find_better = false;
                for (int i = 0; i != pnew.size(); ++i)
                {
                    auto any_dominates = false;
                    for (const auto &z : p)
                    {
                        if (z.strictly_dominates(pnew[i]))
                        {
                            any_dominates = true;
                            if (!find_better)
                            {
                                if (distances.find(pmNew[i]) == distances.end())
                                {
                                    distances[pmNew[i]] = abs(z.distance(pnew[i]));
                                }
                                else
                                {
                                    distances[pmNew[i]] = abs(z.distance(pnew[i])) > abs(distances[pmNew[i]])
                                        ? distances[pmNew[i]]
                                        : abs(z.distance(pnew[i]));
                                }
                            }
                            else
                            {
                                break;
                            }
                        }
                    }

                    if (!any_dominates)
                    {
                        if (success_time.find(pmNew[i]) == success_time.end())
                        {
                            success_time[pmNew[i]] = 1;
                        }
                        else
                        {
                            success_time[pmNew[i]]++;
                        }
                        find_better = true;
                    }
                }


                double best_key;
                if (find_better)
                {
                    int max_success = 0;
                    for (auto const &[key, val] : success_time)
                    {
                        if (val > max_success)
                        {
                            max_success = val;
                            best_key = key;
                        }
                    }
                }
                else
                {
                    double min_distance = std::numeric_limits<double>::max();
                    
                    for (auto const &[key, val] : distances)
                    {
                        if (val < min_distance)
                        {
                            min_distance = val;
                            best_key = key;
                        }
                    }
                }
                if (algorithm_name == "tworate")
                {
                    double best_pm;
                    if (ioh::common::random::real() < 0.5)
                    {
                        best_pm = best_key;
                    }
                    else
                    {
                        if (ioh::common::random::real() < 0.5)
                        {
                            best_pm = pmNew[0];
                        }
                        if (ioh::common::random::real() < 0.5)
                        {
                            best_pm = pmNew[lambda - 1];
                        }
                    }
                    if (best_pm >= 0.5) best_pm = 0.5;
                    if (best_pm < (0.25 / problem->meta_data().n_variables)) best_pm = (0.25 / problem->meta_data().n_variables);
                    r = best_pm * problem->meta_data().n_variables;
                    }
                else if (algorithm_name == "varctrl")
                {
                    r = best_key;
                    if (r >= 0.5) r = 0.5;
                    if (r < (0.25 / problem->meta_data().n_variables)) r = (0.25 / problem->meta_data().n_variables);
                
                    if(find_better) ++c;
                    else c = 1;
                }
                else if (algorithm_name == "logNormal")
                {
                    pm_star = best_key;
                    if (pm_star >= 0.5) pm_star = 0.5;
                    if (pm_star < (0.25 / problem->meta_data().n_variables)) pm_star = (0.25 / problem->meta_data().n_variables);
                }
            }

            /***
             * Update population
             */
            for (const auto &np : pnew)
            {
                auto any_dominates = false;
                std::vector<GSolution> updatep{np};

                for (const auto &z : p)
                {
                    if (z.strictly_dominates(np))
                    {
                        any_dominates = true;
                        break;
                    }

                    if (!z.dominated_by(np))
                        updatep.push_back(z);
                }
                if (!any_dominates)
                    p = updatep;
            }
        }
        return p;
    }


    GSolution eval(const std::vector<int> &x, const std::shared_ptr<ioh::problem::Integer> &problem)
    {
        GSolution s;
        s.x = x;
        s.y[0] = (*problem)(x);
        //std::cout << s.y[0];
        for (size_t i = 1; i < GSolution::S; i++)
        {
            s.y[i] = problem->constraints()[i - 1]->violation();
         //   std::cout << ',' << s.y[i];
        }
        //std::cout << std::endl;
        return s;
    }
};

/***
 * ./main problem_id dimension static/tworate/varctrl/lognormal lambda p budget runs
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
    int budget = std::atoi(argv[6]);
    int runs = std::atoi(argv[7]);

    // int problem_id = 1;
    // int dimension = 100;
    // std::string algorithm_name = "logNormal";
    // int lambda = 2;
    // double pm = 1 / static_cast<double>(dimension);
    // int budget = 100;
    // int runs = 3;

    using namespace ioh::common;
    using namespace ioh::problem;
    random::seed(10);

    bool force_flip = true;


    auto problem = get_problem(problem_id, dimension);
    std::string exp_name = algorithm_name + "L" + argv[4] + "P" + argv[5];
    // std::string exp_name = "test";
    auto logger =
        ioh::logger::Analyzer({ioh::trigger::always}, {ioh::watch::violation}, "/data/yef/MO", exp_name, exp_name, exp_name);


    problem->attach_logger(logger);

    for (int i = 0; i < runs; ++i)
    {
        GSEMO<OptimizationType::MAX, OptimizationType::MAX> opt(budget, force_flip, pm, lambda, algorithm_name);
        auto p = opt(problem);
        //print(p);
        problem->reset();
    }
}
