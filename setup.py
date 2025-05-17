from setuptools import setup, Extension
import pybind11
from setuptools.command.build_ext import build_ext
import os

class BuildExtInSource(build_ext):
    def build_extension(self, ext):
        ext._full_name = ext.name
        ext_path = self.get_ext_fullpath(ext.name)
        # Place the .so in src/pdb_cpp/
        ext_path = os.path.join("src", "pdb_cpp", os.path.basename(ext_path))
        self._set_output(ext, ext_path)
        super().build_extension(ext)

ext_modules = [
    Extension(
        "pdb_cpp",
        ["src/pdb_cpp/model_pybind.cpp", "src/pdb_cpp/Model.cpp", "src/pdb_cpp/format/pdb.cpp"],
        include_dirs=[
            pybind11.get_include(),
            "."
        ],
        language="c++"
    ),
]

setup(
    name="pdb_cpp",
    ext_modules=ext_modules,
    cmdclass={'build_ext': build_ext},
    options={
        'build_ext': {
            'inplace': True,
            'build_lib': 'src/pdb_cpp'
        }
    },
)