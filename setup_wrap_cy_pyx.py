
# >> python setup_wrap_cy_pyx.py build_ext --inplace

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext as _build_ext
from Cython.Build import cythonize
import numpy

include_dirs=[numpy.get_include()]

compile_args = ['/openmp'] #['/O0'] 
link_args = []

extension_kwargs = dict(
    include_dirs=include_dirs,
    extra_compile_args=compile_args,
    extra_link_args=link_args,
    language='c',
    )


class build_ext(_build_ext):
    def build_extensions(self):
        super().build_extensions()

    def initialize_options(self):
        super().initialize_options()
        if self.distribution.ext_modules == None:
            self.distribution.ext_modules = []

        extensions=[
            Extension("*", sources=["./myfempy/core/solver/assemblersymm_cython_v4.pyx"], **extension_kwargs),
            Extension("*", sources=["./myfempy/core/solver/assemblerfull_cython_v4.pyx"], **extension_kwargs),
            Extension("*", sources=["./myfempy/core/shapes/line2_tasks.pyx"], **extension_kwargs),
            Extension("*", sources=["./myfempy/core/shapes/line3_tasks.pyx"], **extension_kwargs),
            Extension("*", sources=["./myfempy/core/shapes/tria3_tasks.pyx"], **extension_kwargs),
            Extension("*", sources=["./myfempy/core/shapes/tria6_tasks.pyx"], **extension_kwargs),
            Extension("*", sources=["./myfempy/core/shapes/quad4_tasks.pyx"], **extension_kwargs),
            Extension("*", sources=["./myfempy/core/shapes/quad8_tasks.pyx"], **extension_kwargs),
            Extension("*", sources=["./myfempy/core/shapes/tetr4_tasks.pyx"], **extension_kwargs),
        ]
        self.distribution.ext_modules.extend(cythonize(extensions))