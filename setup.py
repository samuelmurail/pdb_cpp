import sys

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import pybind11


def _cpp_args():
    if sys.platform == "win32":
        return [
            "/O2",
            "/fp:fast",
            "/std:c++17",
            "/wd4244",
            "/wd4267",
            "/wd4305",
            "/wd4996",
        ]
    return ["-O3", "-ffast-math", "-std=c++17"]

ext_modules = [
    Extension(
        "pdb_cpp.core",
        ["src/pdb_cpp/_core/pybind.cpp",
         "src/pdb_cpp/_core/Coor.cpp",
         "src/pdb_cpp/_core/Model.cpp",
         "src/pdb_cpp/_core/format/pdb.cpp",
         "src/pdb_cpp/_core/format/encode.cpp",
         "src/pdb_cpp/_core/format/mmcif.cpp",
         "src/pdb_cpp/_core/select.cpp",
         "src/pdb_cpp/_core/sequence.cpp",
         "src/pdb_cpp/_core/align.cpp",
         "src/pdb_cpp/_core/seq_align.cpp",
         "src/pdb_cpp/_core/data/residue.cpp",
         "src/pdb_cpp/_core/TMalign_wrapper.cpp"],
        include_dirs=[
            pybind11.get_include(),
            "src/pdb_cpp/_core",
        ],
        extra_compile_args=_cpp_args(),
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