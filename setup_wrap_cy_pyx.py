
# >> python setup_cython_wrap.py build_ext --inplace

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext as _build_ext
from Cython.Build import cythonize
import numpy

include_dirs=[numpy.get_include()]

compile_args = ['/O0']
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
            # Extension("*", sources=["./myfempy/expe/utilities.pyx"], **extension_kwargs),
            Extension("*", sources=["./myfempy/expe/asmb_cython/assembler_cython.pyx"], **extension_kwargs),
        ]
        self.distribution.ext_modules.extend(cythonize(extensions))