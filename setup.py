import os
from setuptools import setup

from pybind11.setup_helpers import Pybind11Extension, build_ext

__version__ = "0.0.1"
os.environ["CC"] = "g++"
os.environ["CXX"] = "g++"

ext = Pybind11Extension(
    "gsemo.gsemocpp",
    ["src/problems.cpp", "src/interface.cpp"],
    include_dirs=[
        "include",
        "src",
        "external/fmt/include",
        "external/cxxopts/include",
        "external/json/include",
    ],
    cxx_std=17,
)

ext._add_cflags(["-O3"])

setup(
    name="GSEMO",
    version=__version__,
    author="Jacob de Nobel, Furong Ye",
    long_description="",
    ext_modules=[ext],
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
)
