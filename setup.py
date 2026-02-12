import sys

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import pybind11


def _cpp_args():
    if sys.platform == "win32":
        return ["/std:c++17"]
    return ["-std=c++17"]

ext_modules = [
    Extension(
        "pdb_cpp.core",
        ["src/pdb_cpp/pybind.cpp",
         "src/pdb_cpp/Coor.cpp",
         "src/pdb_cpp/Model.cpp",
         "src/pdb_cpp/format/pdb.cpp",
         "src/pdb_cpp/format/encode.cpp",
         "src/pdb_cpp/format/mmcif.cpp",
         "src/pdb_cpp/select.cpp",
         "src/pdb_cpp/sequence.cpp",
         "src/pdb_cpp/align.cpp",
         "src/pdb_cpp/seq_align.cpp",
         "src/pdb_cpp/data/residue.cpp",
         "src/pdb_cpp/TMalign_wrapper.cpp"],
        include_dirs=[
            pybind11.get_include(),
            "src/pdb_cpp",
        ],
        extra_compile_args=["-O3", "-ffast-math"] + _cpp_args(),
        language="c++"
    ),
]

setup(
    name="pdb_cpp",
    package_dir={"": "src"},
    packages=["pdb_cpp"],
    ext_modules=ext_modules,
    cmdclass={'build_ext': build_ext},
)