#include "gsemo.hpp"
#include "pybind11/numpy.h"
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(gsemocpp, m)
{
    using namespace ioh::common;

    m.def("get_problem", &get_problem);
    m.def("print_front", &print<OptimizationType::MAX, OptimizationType::MAX>);

    using Solution = MultiSolution<OptimizationType::MAX, OptimizationType::MAX>;
    py::class_<Solution>(m, "MultiSolution")
        .def_readonly("x", &Solution::x)
        .def_readonly("y", &Solution::y)
        .def_readwrite("pm", &Solution::y)
        .def_readwrite("l", &Solution::y)
        .def("__eq__", &Solution::operator==)
        .def("dominated_by", &Solution::dominated_by)
        .def("strictly_dominates", &Solution::strictly_dominates)
        .def("distance", &Solution::distance)
        .def("print", &Solution::print)
        .def("eval", &Solution::eval)
        .def("__repr__", [](const Solution& self) {return fmt::format("({})", fmt::join(self.y, ", "));})
        ;

    using GESMO_ = GSEMO<OptimizationType::MAX, OptimizationType::MAX>;
    py::class_<GESMO_>(m, "GSEMO")
        .def(
            py::init<int, bool, double, int, std::string, int>(),
            py::arg("budget") = 6000,
            py::arg("force_flip") = true,
            py::arg("pm") = 0.01,
            py::arg("lambda") = 1,
            py::arg("algorithm_name") = "static",
            py::arg("verbose_rate") = 0)
        .def("__call__", &GESMO_::operator())
        .def("evaluate_problem", [](GESMO_ &self, const int id, const int dim, const bool logged, const int nreps, const std::string& path)
             {
            auto problem = get_problem(id, dim);
            const std::string exp_name = fmt::format("{}L{}PM0{}", self.algorithm_name, self.strategy->lambda, self.strategy->pm0);
            auto logger =
                ioh::logger::Analyzer({ioh::trigger::always}, {ioh::watch::violation}, path, exp_name, exp_name, exp_name);   
            
            if (logged){
                problem->attach_logger(logger);         
            }
            std::vector<std::vector<GESMO_::GSolution>> res(nreps);
            for (int i = 0; i < nreps; ++i){
                res[i] = self(problem);
                problem->reset();
            }
            return res; });
}