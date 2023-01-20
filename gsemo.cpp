#include "problems_gecco23.hpp"

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

    // @Furong: Why is this templated?
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

    //! Sampling n different indexes from length man
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
        {
            /// If n is smaller than m/2, we sample indexes repeatly until getting different values.
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
    double pm;
    int l;

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
            [&]
            {
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
            [&]
            {
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
            [&]
            {
                dis += operators::distance_<Args>(y[i], other.y[i]);
                i++;
            }(),
            ...);
        return dis;
    }

    void distance_to_front(const std::vector<MultiSolution<Args...>> &other)
    {
    }

    void print() const
    {
        std::cout << "(";
        for (const auto &fi : y)
            std::cout << fi << ",";
        std::cout << ")";
    }

    void eval(const std::shared_ptr<ioh::problem::IntegerSingleObjective> &problem)
    {
        y[0] = (*problem)(x);
        for (size_t i = 1; i < S; i++)
            y[i] = problem->constraints()[i - 1]->violation();

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
    const auto n = std::max(m, d(ioh::common::random::GENERATOR));
    for (const auto &pos : ioh::common::random::integers(n, 0, (int)x.size() - 1))
        x[pos] = abs(x[pos] - 1);
}

namespace adaptation
{
    struct PopulationStats
    {
        bool has_non_dominated_solution;
        std::vector<bool> successes;
        std::vector<double> distances;

        PopulationStats(size_t n) : has_non_dominated_solution(false), successes(n, false), distances(n, 0.) {}

        template <typename SolutionType>
        PopulationStats(const std::vector<SolutionType> &pareto_front, const std::vector<SolutionType> &new_population) : PopulationStats(new_population.size())
        {
            for (size_t i = 0; i < new_population.size(); i++)
            {
                auto any_dominates = false;
                for (const auto &z : pareto_front)
                {
                    distances[i] += abs(z.distance(new_population[i]));
                    any_dominates |= z.strictly_dominates(new_population[i]);
                }
                successes[i] = !any_dominates;
                has_non_dominated_solution |= !any_dominates;
            }
        }
    };

    template <typename SolutionType>
    struct Strategy
    {
        double pm0;
        double pm;
        int lambda;
        int n;

        Strategy(const double pm, const int lambda) : pm0(pm), pm(pm), lambda(lambda) {}

        virtual double generate_pm(const int) const
        {
            return pm0;
        };

        virtual int generate_l(const double pm) const
        {
            std::binomial_distribution<> d{n, pm};
            int l = d(ioh::common::random::GENERATOR);
            while (l < 1)
                l = d(ioh::common::random::GENERATOR);
            return l;
        }

        virtual void setup_problem(const std::shared_ptr<ioh::problem::IntegerSingleObjective> &problem)
        {
            n = problem->meta_data().n_variables;
        }

        virtual void adapt(const std::vector<SolutionType> &pareto_front, const std::vector<SolutionType> &new_population) {}
    };

    template <typename SolutionType>
    using Static = Strategy<SolutionType>;

    template <typename SolutionType>
    struct TwoRate : Strategy<SolutionType>
    {
        using Strategy<SolutionType>::Strategy;

        double generate_pm(const int i) const override
        {
            if (i < this->lambda / 2)
                return r / 2.0 / this->n;
            return r * 2.0 / this->n;
        };

        void setup_problem(const std::shared_ptr<ioh::problem::IntegerSingleObjective> &problem) override
        {
            Strategy<SolutionType>::setup_problem(problem);
            r = this->pm0 * this->n;
        }

        virtual void adapt(const std::vector<SolutionType> &pareto_front, const std::vector<SolutionType> &new_population) override
        {
            // TODO: this is was not defined according to the paper, I changed it.
            const auto &stats = PopulationStats(pareto_front, new_population);
            const size_t half = stats.successes.size() / 2;
            double s;
            if (stats.has_non_dominated_solution)
            {
                const auto s1 = std::accumulate(stats.successes.begin(), stats.successes.begin() + half, 0);
                const auto s2 = std::accumulate(stats.successes.begin() + half, stats.successes.end(), 0);
                s = s1 > s2 ? 0.75 : 0.25;
            }
            else
            {
                const auto min_element = std::distance(stats.distances.begin(), std::min_element(stats.distances.begin(), stats.distances.end()));
                s = min_element < half ? 0.75 : 0.25;
            }
            if (ioh::common::random::real() <= s) {
                r = std::max(r / 2., 0.5);
            } else {
                r = std::min(2. * r, this->n / 4.);
            }
        }

    protected:
        double r;
    };

    template <typename SolutionType>
    struct LogNormal : Strategy<SolutionType>
    {
        using Strategy<SolutionType>::Strategy;

        double generate_pm(const int i) const override
        {
            static std::normal_distribution<> standard_norm{0, 1};

            const auto r = standard_norm(ioh::common::random::GENERATOR);

            return 1.0 / (1.0 + (((1.0 - this->pm) / this->pm) * exp(0.22 * r)));
        }
        virtual void adapt(const std::vector<SolutionType> &pareto_front, const std::vector<SolutionType> &new_population) override
        {
            auto stats = PopulationStats(pareto_front, new_population);
            size_t i;
            if (!stats.has_non_dominated_solution){
                i = std::distance(stats.distances.begin(), std::min_element(stats.distances.begin(), stats.distances.end()));
            } else {
                i = std::distance(stats.successes.begin(), std::find(stats.successes.begin(), stats.successes.end(), true));
            }
            this->pm = std::max(0.25 / this->n, std::min(0.5, new_population[i].pm));
        }
    };

    template <typename SolutionType>
    struct VarCntrl : TwoRate<SolutionType>
    {
        using TwoRate<SolutionType>::TwoRate;

        int generate_l(const double pm) const
        {
            const double dn = static_cast<double>(this->n);
            std::normal_distribution<> d{this->r, this->r * (1 - this->r / dn) * pow(F, c)};

            double l = d(ioh::common::random::GENERATOR);
            while (l < 1 || l > this->n / 2)
                l = d(ioh::common::random::GENERATOR);

            return static_cast<double>(std::round(l));
        }

        virtual void adapt(const std::vector<SolutionType> &pareto_front, const std::vector<SolutionType> &new_population) override
        {
            auto stats = PopulationStats(pareto_front, new_population);
            size_t i;
            if (!stats.has_non_dominated_solution){
                c = 0;
                i = std::distance(stats.distances.begin(), std::min_element(stats.distances.begin(), stats.distances.end()));
            } else {
                c++;
                i = std::distance(stats.successes.begin(), std::find(stats.successes.begin(), stats.successes.end(), true));
            }

            this->r = new_population[i].l;
        }

    private:
        double F = 0.98;
        int c = 0;
    };

    template <typename T>
    std::unique_ptr<adaptation::Strategy<T>> get(const std::string &name, const double pm, const int lambda)
    {
        if (name == "static")
            return std::make_unique<adaptation::Static<T>>(pm, lambda);
        else if (name == "TwoRate")
            return std::make_unique<adaptation::TwoRate<T>>(pm, lambda);
        else if (name == "varctrl")
            return std::make_unique<adaptation::VarCntrl<T>>(pm, lambda);
        else if (name == "logNormal")
            return std::make_unique<adaptation::LogNormal<T>>(pm, lambda);
        throw std::invalid_argument(name + " is an unknown adaptation strategy type");
    }
}

template <ioh::common::OptimizationType... Args>
struct GSEMO
{
    using GSolution = MultiSolution<Args...>;

    int budget;
    bool force_flip; // TODO: this is now unused
    int verbose_rate;
    std::string algorithm_name;
    std::unique_ptr<adaptation::Strategy<GSolution>> strategy;

    GSEMO(const int b = 6000,
          const bool ff = true,
          const double pm = 0.01,
          const int lambda = 1,
          const std::string algorithm_name = "static",
          const int v = 0)
        : budget(b),
          force_flip(ff),
          verbose_rate(v),
          algorithm_name(algorithm_name),
          strategy(adaptation::get<GSolution>(algorithm_name, pm, lambda))
    {
    }

    std::vector<GSolution> operator()(const std::shared_ptr<ioh::problem::IntegerSingleObjective> &problem)
    {
        std::vector<GSolution> pareto_front{1};
        pareto_front[0].x = ioh::common::random::integers(problem->meta_data().n_variables, 0, 1);
        pareto_front[0].eval(problem);
        --budget;

        strategy->setup_problem(problem);

        while (budget > 0)
        {
            auto new_population = std::vector<GSolution>(strategy->lambda);

            for (int i = 0; i < strategy->lambda; i++)
            {
                auto& candidate = new_population[i];

                // Select new candidate
                const auto ri = ioh::common::random::integer(0, pareto_front.size() - 1);
                candidate.x = pareto_front[ri].x;
                candidate.pm = strategy->generate_pm(i);
                candidate.l = strategy->generate_l(candidate.pm);

                // Mutate
                bitflip_index(candidate.x, candidate.l);

                // Evaluate
                candidate.eval(problem);

                if (--budget <= 0 || problem->state().optimum_found)
                    break;

                if (verbose_rate and i % verbose_rate == 0)
                    print(pareto_front);
            }

            strategy->adapt(pareto_front, new_population);

            for (const auto &np : new_population)
            {
                auto any_dominates = false;
                std::vector<GSolution> non_dominated_solutions{np};

                for (const auto &z : pareto_front)
                {
                    if (z.strictly_dominates(np))
                    {
                        any_dominates = true;
                        break;
                    }

                    if (!z.dominated_by(np))
                        non_dominated_solutions.push_back(z);
                }
                if (!any_dominates)
                    pareto_front = non_dominated_solutions;
            }
        }
        return pareto_front;
    }
};

/***
 * ./main problem_id dimension static/tworate/varctrl/lognormal lambda p budget runs
 */
int main(int argc, char *argv[])
{
    // if (argc < 6)
    // {
    //     std::cerr << "Some parameters are missing" << std::endl;
    // }

    // int problem_id = std::atoi(argv[1]);
    // int dimension = std::atoi(argv[2]);
    // std::string algorithm_name = argv[3];
    // int lambda = std::atoi(argv[4]);
    // double pm = std::atoi(argv[5]) / static_cast<double>(dimension);
    // int budget = std::atoi(argv[6]);
    // int runs = std::atoi(argv[7]);

    int problem_id = 1;
    int dimension = 100;
    int lambda = 2;
    double pm = 1 / static_cast<double>(dimension);
    int budget = 100000;
    int runs = 1;

    using namespace ioh::common;
    random::seed(10);

    bool force_flip = true;

    auto problem = get_problem(problem_id, dimension);
    
    std::string algorithm_name = "static";
    std::string exp_name = algorithm_name + "L" + argv[4] + "P" + argv[5];
    // std::string exp_name = "test";
    // auto logger =
    //     ioh::logger::Analyzer({ioh::trigger::always}, {ioh::watch::violation}, "/data/yef/MO", exp_name, exp_name, exp_name);

    // problem->attach_logger(logger);

    for (int i = 0; i < runs; ++i)
    {
        GSEMO<OptimizationType::MAX, OptimizationType::MAX> opt(budget, force_flip, pm, lambda, algorithm_name);
        auto p = opt(problem);
        print(p);
        std::cout << p.size() << std::endl;
        std::cout << *problem << std::endl;
        problem->reset();
    }
}
