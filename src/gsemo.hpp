#include "problems.hpp"

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
                dis += pow(operators::distance_<Args>(y[i], other.y[i]),2);
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
        for (size_t i = 1; i < S; i++) {
            y[i] = problem->constraints()[i - 1]->violation();
        }
    }
};

template <ioh::common::OptimizationType... Args>
void print(const std::vector<MultiSolution<Args...> > &p)
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
        // Make sure that the references point are ranked by the first obj.
        std::vector<SolutionType> pareto_reference;
        SolutionType hv_reference;
        int adapt_metric; // 1: hypervolume 2: IGD 3: NOPareto

        Strategy(const double pm, const int lambda, const int adapt_metric = 1) : pm0(pm), pm(pm), lambda(lambda), adapt_metric(adapt_metric) {}

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

        virtual void setup_problem(const std::shared_ptr<ioh::problem::IntegerSingleObjective> &problem, const std::vector<SolutionType> &pareto_front, const SolutionType & hv)
        {
            this->n = problem->meta_data().n_variables;
            this->pareto_reference = pareto_front;
            this->hv_reference = hv;
            // std::sort(pareto_reference.begin(), pareto_reference.end(), [](const SolutionType& a, const SolutionType& b) {return a.y[0] < b.y[0];});
        }

        virtual void adapt(const std::vector<SolutionType> &pareto_front, const std::vector<SolutionType> &new_population) {}

        double pareto_metric(const std::vector<SolutionType> &new_population) {
            auto tmp_pareto_front = new_population;
            std::sort(tmp_pareto_front.begin(),tmp_pareto_front.end());
            switch (adapt_metric)
            {
                case 1:
                    return hypervolume(tmp_pareto_front);
                    break;
                case 2:
                    return -inverted_generational_distance(tmp_pareto_front);
                    break;
                case 3:
                    return number_of_pareto(tmp_pareto_front);
                    break;
                default:
                    return -1;
                    break;
            }   
        }

        double hypervolume(const std::vector<SolutionType> &new_population) 
        {
            // std::sort(new_population.begin(), new_population.end());
            double volume = 0;
            double last_x = this->hv_reference.y[0];
            for (const auto & p : new_population)
            {   
                volume += (p.y[1] - this->hv_reference.y[1]) * (p.y[0] - last_x);
                last_x = p.y[0];
            }
            return volume;
        }

        double inverted_generational_distance(const std::vector<SolutionType> &new_population) 
        {
            // std::sort(new_population.begin(), new_population.end(), [](const SolutionType& a, const SolutionType& b) {return a.y[0] < b.y[0];});
            double igd = 0;
            for (size_t i = 0; i != this->pareto_reference.size(); ++i)
            {
                double tmp = n * n  + 1;
                for(size_t j = 0; j != new_population.size(); ++j) 
                {
                    if ( pow((new_population[j].y[0] - pareto_reference[i].y[0]),2) > tmp  ) { break; } 
                    double d = pareto_reference[i].euclidean_distance(new_population[j]);
                    if (tmp < d) { tmp = d;}
                }
                igd += pow(tmp,2);
            }
            return igd / this->pareto_reference.size();
        }

        double number_of_pareto(const std::vector<SolutionType> &new_population) 
        {
            // std::sort(new_population.begin(), new_population.end());
            double found = 0;
            size_t j = 0;
            for (size_t i = 0; i != this->pareto_reference.size();)
            {
                if(pareto_reference[i].y[0] == new_population[j].y[0]) {
                    if(pareto_reference[i] == new_population[j]) {found += 1.0;}
                    ++i;
                    ++j;
                }
                else if(pareto_reference[i].y[0] < new_population[j].y[0]) {
                    ++i;
                } 
                else if(pareto_reference[i].y[0] > new_population[j].y[0]) {
                    ++j;
                }
                if (j == new_population.size()) break;
            }
            return found;
        }

        
        std::vector<double> population_stats(const std::vector<SolutionType> &pareto_front, const std::vector<SolutionType> &new_population)
        {
            
            std::vector<double> metrics(new_population.size(), -1.0);
            auto tmp_pareto_front = pareto_front;
            for (size_t i = 0; i < new_population.size(); i++)
            {
                tmp_pareto_front.push_back(new_population[i]);
                metrics[i] = this->pareto_metric(tmp_pareto_front);
                tmp_pareto_front.pop_back();
            }
            return metrics;
        }
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

        void setup_problem(const std::shared_ptr<ioh::problem::IntegerSingleObjective> &problem, const std::vector<SolutionType> &pareto_front, const SolutionType & hv) override
        {
            Strategy<SolutionType>::setup_problem(problem,pareto_front,hv);
            r = this->pm0 * this->n;
        }

        virtual void adapt(const std::vector<SolutionType> &pareto_front, const std::vector<SolutionType> &new_population) override
        {
            // TODO: this is was not defined according to the paper, I changed it.
            auto metrics = this->population_stats(pareto_front, new_population);
            const size_t half = new_population.size() / 2;
            double s;
            size_t i = std::distance(metrics.begin(), std::max_element(metrics.begin(), metrics.end()));
            s = i < half ? 0.75 : 0.25;
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
            auto metrics = this->population_stats(pareto_front, new_population);
            size_t i;
            i = std::distance(metrics.begin(), std::max_element(metrics.begin(), metrics.end()));
            this->pm = std::max(0.25 / this->n, std::min(0.5, new_population[i].pm));
        }
    };

    template <typename SolutionType>
    struct VarCntrl : TwoRate<SolutionType>
    {
        using TwoRate<SolutionType>::TwoRate;

        int generate_l(const double pm) const override
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
            auto metrics = this->population_stats(pareto_front, new_population);
            size_t i;
            i = std::distance(metrics.begin(), std::max_element(metrics.begin(), metrics.end()));
            if (this->r == new_population[i].l) 
            {
                c += 1;
            }
            else 
            {
                c = 0;
                this->r = new_population[i].l;
            }
        }

        void setup_problem(const std::shared_ptr<ioh::problem::IntegerSingleObjective> &problem, const std::vector<SolutionType> &pareto_front, const SolutionType & hv) override
        {
            TwoRate<SolutionType>::setup_problem(problem,pareto_front,hv);
            c = 0;
        }

    private:
        double F = 0.98;
        int c = 0;
    };
    
    template <typename SolutionType>
    struct NormEA : TwoRate<SolutionType>
    {
        using TwoRate<SolutionType>::TwoRate;

        int generate_l(const double pm) const override
        {
            const double dn = static_cast<double>(this->n);
            std::normal_distribution<> d{this->r, this->r * (1 - this->r / dn)};

            double l = d(ioh::common::random::GENERATOR);
            while (l < 1 || l > this->n / 2)
                l = d(ioh::common::random::GENERATOR);

            return static_cast<double>(std::round(l));
        }

        virtual void adapt(const std::vector<SolutionType> &pareto_front, const std::vector<SolutionType> &new_population) override
        {
            auto metrics = this->population_stats(pareto_front, new_population);
            size_t i;
            i = std::distance(metrics.begin(), std::max_element(metrics.begin(), metrics.end()));
            this->r = new_population[i].l;
        }

        void setup_problem(const std::shared_ptr<ioh::problem::IntegerSingleObjective> &problem, const std::vector<SolutionType> &pareto_front, const SolutionType & hv) override
        {
            TwoRate<SolutionType>::setup_problem(problem,pareto_front,hv);
        }
    };


    template <typename T>
    std::unique_ptr<adaptation::Strategy<T> > get(const std::string &name, const double pm, const int lambda, const int adapt_metric)
    {
        if (name == "static")
            return std::make_unique<adaptation::Static<T> >(pm, lambda);
        else if (name == "TwoRate")
            return std::make_unique<adaptation::TwoRate<T> >(pm, lambda, adapt_metric);
        else if (name == "varctrl")
            return std::make_unique<adaptation::VarCntrl<T> >(pm, lambda, adapt_metric);
        else if (name == "normea")
            return std::make_unique<adaptation::NormEA<T> >(pm, lambda, adapt_metric);
        else if (name == "logNormal")
            return std::make_unique<adaptation::LogNormal<T> >(pm, lambda, adapt_metric);
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
    std::unique_ptr<adaptation::Strategy<GSolution> > strategy;

    GSEMO(const int b = 6000,
          const bool ff = true,
          const double pm = 0.01,
          const int lambda = 1,
          const int adapt_metric = 1,
          const std::string algorithm_name = "static",
          const int v = 0)
        : budget(b),
          force_flip(ff),
          verbose_rate(v),
          algorithm_name(algorithm_name),
          strategy(adaptation::get<GSolution>(algorithm_name, pm, lambda, adapt_metric))
    {
    }

    static bool compare(const GSolution & a, const GSolution &b) {
        return a.y[0] < b.y[0];
    }

    std::vector<GSolution> operator()(const std::shared_ptr<ioh::problem::IntegerSingleObjective> &problem)
    {
        // Preparing pareto_front_star and hv_reference
        auto pareto_star = std::vector<GSolution>(pre_pareto_y.size());
        for (size_t i = 0; i != pre_pareto_y.size(); ++i)
        {
            pareto_star[i].y[0] = pre_pareto_y[i][0];
            pareto_star[i].y[1] = pre_pareto_y[i][1];
        }
        GSolution hv;
        hv.y[0] = pre_hv_reference[0];
        hv.y[1] = pre_hv_reference[1]; 

        std::vector<GSolution> pareto_front{1};
        pareto_front[0].x = ioh::common::random::integers(problem->meta_data().n_variables, 0, 1);
        pareto_front[0].eval(problem);
        strategy->setup_problem(problem,pareto_star,hv);

        while (problem->state().evaluations < budget)
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
                bitflip(candidate.x, candidate.l);

                // Evaluate
                candidate.eval(problem);

                if (problem->state().evaluations > budget)
                    break;

                if (verbose_rate and i % verbose_rate == 0)
                    print(pareto_front);
            }
            
            // Sort the current 
            // std::sort(new_population.begin(), new_population.end())
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
            if (pareto_front.size() ==  pareto_star.size()) 
            {
                bool hit_all_pareto = true;
                std::sort(pareto_front.begin(), pareto_front.end(),compare);
                for (size_t i = 0; i != pareto_star.size(); ++i)
                {
                    if ( (pareto_front[i].y[0] != pareto_star[i].y[0] ) || (pareto_front[i].y[1] != pareto_star[i].y[1]) )
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
            
        }
        return pareto_front;
    }
};
