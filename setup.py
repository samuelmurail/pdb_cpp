from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import pybind11

ext_modules = [
    Extension(
        "pdb_cpp.pdb_cpp",
        ["src/pdb_cpp/model_pybind.cpp", "src/pdb_cpp/Model.cpp", "src/pdb_cpp/format/pdb.cpp"],
        include_dirs=[
            pybind11.get_include(),
            "src/pdb_cpp"
        ],
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