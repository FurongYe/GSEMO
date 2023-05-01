#include "problems.hpp"
#include "math.h"

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
    std::array<double, S> p_y;
    std::vector<int> x;

    double previous_dis = -1;
    double dis = -1;
    int c = 0;
    double pm = 0.01;
    double r = 1;
    int l;

    bool operator==(const MultiSolution<Args...> &other) const
    {
        for (size_t i = 0; i < S; i++)
            if (y[i] != other.y[i])
                return false;
        return true;
    }

    bool operator<(const MultiSolution<Args...> &other) const
    {
        return (y[0] < other.y[0]);
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

    double euclidean_distance(const MultiSolution<Args...> &other) const
    {
        double dis;
        int i = 0;
        (
            [&]
            {
                dis += pow(operators::distance_<Args>(y[i], other.y[i]), 2);
                i++;
            }(),
            ...);
        dis = sqrt(dis);
        return dis;
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
        {
            y[i] = problem->constraints()[i - 1]->violation();
        }
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

void bitflip_index(std::vector<int> &x, const int flip_bits_number)
{
    std::vector<size_t> flip_index;
    operators::sampleNFromM(flip_index, flip_bits_number, x.size());
    for (const auto &pos : flip_index)
        x[pos] = abs(x[pos] - 1);
}

void bitflip(std::vector<int> &x, const int n)
{
    for (const auto &pos : ioh::common::random::integers(n, 0, (int)x.size() - 1))
        x[pos] = abs(x[pos] - 1);
}

namespace adaptation
{

    template <typename SolutionType>
    struct Strategy
    {
        double pm0;
        double pm;
        int lambda;
        int n;
        int generation = 1;

        Strategy(const double pm, const int lambda) : pm0(pm), pm(pm), lambda(lambda) {}

        virtual double generate_pm(const int i, const std::vector<SolutionType> &a) const
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

        virtual int generate_l_n(const double r, const int c) const
        {
           return 1;
        }

        virtual void setup_problem(const std::shared_ptr<ioh::problem::IntegerSingleObjective> &problem)
        {
            this->n = problem->meta_data().n_variables;
        }

        virtual void adapt(const std::vector<SolutionType> &pareto_front,  std::vector<SolutionType> &new_population) {}
    };

    template <typename SolutionType>
    struct AGSEMO_Strategy : Strategy<SolutionType>
    {
        using Strategy<SolutionType>::Strategy;
        
        double F = 0.98;

        // parameters for the standard bit mutation, to be applied to the other x.
        double r = 1;

        double generate_pm(const int i, const std::vector<SolutionType> &pareto_front) const override
        {
            
            // for the other x
            double tmp = pareto_front[i].r /(double)this->n;
            static std::normal_distribution<> standard_norm{0, 1};

            const auto r = standard_norm(ioh::common::random::GENERATOR);

            tmp =  1.0 / (1.0 + (((1.0 - tmp) / tmp) * exp(0.22 * r)));


            tmp = tmp > 0.5 ? 0.5 : tmp;
            tmp = tmp < 1.0 / this->n ? 1.0 / this->n : tmp;
            return tmp;
        }

        int generate_l(const double pm) const override
        {
            
                std::binomial_distribution<> d{this->n, pm};
                int l = d(ioh::common::random::GENERATOR);
                while (l < 1)
                    l = d(ioh::common::random::GENERATOR);
                return l;
            
        }

        int generate_l_n(const double r, const int c) const override {
            
            std::normal_distribution<> d{r, std::sqrt(r * (1 - r / this->n) * pow(F, c))};

            double l = std::round(d(ioh::common::random::GENERATOR));
            while (l < 1 || l > this->n / 2)
                l = d(ioh::common::random::GENERATOR);

            return static_cast<double>(l);
        }

        virtual void adapt(const std::vector<SolutionType> &pareto_front, std::vector<SolutionType> &new_population) override
        {
            for (size_t i = 0; i < new_population.size(); i++)
            {
                
                if (new_population[i].r == new_population[i].l)
                {
                    new_population[i].c += 1.0;
                }
                else{
                    new_population[i].r = new_population[i].l;
                    new_population[i].c = 0;
                }
            }
        }
    };

    template <typename T>
    std::unique_ptr<adaptation::Strategy<T>> getAGSEMO_Strategy(const double pm, const int lambda)
    {
        return std::make_unique<adaptation::AGSEMO_Strategy<T>>(pm, lambda);
    }
}

template <ioh::common::OptimizationType... Args>
struct AGSEMO
{
    using GSolution = MultiSolution<Args...>;
    int budget;
    bool force_flip; // TODO: this is now unused
    int verbose_rate;
    std::string algorithm_name;
    std::unique_ptr<adaptation::Strategy<GSolution>> strategy;
    std::vector<GSolution> pareto_front;

    AGSEMO(const int b = 6000,
           const bool ff = true,
           const double pm = 0.01,
           const int lambda = 1,
           const std::string algorithm_name = "AGSEMO",
           const int v = 0)
        : budget(b),
          force_flip(ff),
          verbose_rate(v),
          algorithm_name(algorithm_name),
          strategy(adaptation::getAGSEMO_Strategy<GSolution>(pm, lambda))
    {
    }

    static bool compare(const GSolution &a, const GSolution &b)
    {
        return a.y[0] < b.y[0];
    }

    void pareto_dis(std::vector<GSolution> & pareto_front) {
        double dis, dis1;
        std::sort(pareto_front.begin(), pareto_front.end(), compare);
        
        dis = std::sqrt(pow(pareto_front[0].y[0] - pareto_front[1].y[0], 2) + pow(pareto_front[0].y[1] - pareto_front[1].y[1], 2));
        pareto_front[0].previous_dis = pareto_front[0].dis;
        pareto_front[0].dis = dis;

        for (size_t i = 1; i < pareto_front.size() - 1; ++i)
        {
            dis = std::sqrt(pow(pareto_front[i - 1].y[0] - pareto_front[i].y[0], 2) + pow(pareto_front[i - 1].y[1] - pareto_front[i].y[1], 2));
            dis1 = std::sqrt(pow(pareto_front[i + 1].y[0] - pareto_front[i].y[0], 2) + pow(pareto_front[i + 1].y[1] - pareto_front[i].y[1], 2));
            dis = std::max(dis, dis1);
            
            pareto_front[i].previous_dis = pareto_front[i].dis;
            pareto_front[i].dis = dis;
        }

        dis = std::sqrt(pow(pareto_front[pareto_front.size() - 2].y[0] - pareto_front[pareto_front.size() - 1].y[0], 2) + pow(pareto_front[pareto_front.size() - 2].y[1] - pareto_front[pareto_front.size() - 1].y[1], 2));
        pareto_front[pareto_front.size() - 1].previous_dis = pareto_front[pareto_front.size() - 1].dis;
        pareto_front[pareto_front.size() - 1].dis = dis;
    }

    std::vector<GSolution> operator()(const std::shared_ptr<ioh::problem::IntegerSingleObjective> &problem)
    {
        // Preparing pareto_front_star
        // std::cout << "run"<< std::endl;
        auto pareto_star = std::vector<GSolution>(pre_pareto_y.size());
        for (size_t i = 0; i != pre_pareto_y.size(); ++i)
        {
            pareto_star[i].y[0] = pre_pareto_y[i][0];
            pareto_star[i].y[1] = pre_pareto_y[i][1];
        }

        pareto_front = std::vector<GSolution>{1};
        pareto_front[0].x = ioh::common::random::integers(problem->meta_data().n_variables, 0, 1);
        pareto_front[0].eval(problem);
        strategy->setup_problem(problem);

        while (problem->state().evaluations < budget)
        {
            auto new_population = std::vector<GSolution>(strategy->lambda);

            for (int i = 0; i < strategy->lambda; i++)
            {
                auto &candidate = new_population[i];
                int ri;
                ri = ioh::common::random::integer(0, pareto_front.size() - 1);
                
                candidate.x = pareto_front[ri].x;
                candidate.p_y = pareto_front[ri].y;
                if ((ri == 0) || (ri == pareto_front.size()-1)) {
                    candidate.l = strategy->generate_l_n(candidate.r,candidate.c);
                    candidate.pm = candidate.l / (double)strategy->n;
                }
                else{
                    candidate.pm = strategy->generate_pm(ri, pareto_front);
                    candidate.l = strategy->generate_l(candidate.pm);
                }
                // std::cout << candidate.l << std::endl;

                // Mutate
                bitflip(candidate.x, candidate.l);

                // Evaluate
                candidate.eval(problem);

                if (problem->state().evaluations > budget)
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

            // Check if all the pareto solutions are found.
            // We know that both lists are sorted.
            if (pareto_front.size() == pareto_star.size())
            {
                bool hit_all_pareto = true;
                std::sort(pareto_front.begin(), pareto_front.end(), compare);
                for (size_t i = 0; i != pareto_star.size(); ++i)
                {
                    if ((pareto_front[i].y[0] != pareto_star[i].y[0]) || (pareto_front[i].y[1] != pareto_star[i].y[1]))
                    {
                        hit_all_pareto = false;
                        break;
                    }
                }
                if (hit_all_pareto)
                {
                    break;
                }
            }

            pareto_dis(pareto_front);
        }
        return pareto_front;
    }
};
