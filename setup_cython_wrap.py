
# >> python setup_cython_wrap.py build_ext --inplace

from setuptools import Extension, setup
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


ext_modules=[
    Extension("*", sources=["./myfempy/experimental/utilities.pyx"], **extension_kwargs),
    Extension("*", sources=["./myfempy/experimental/asmb_cython/assembler_cython.pyx"], **extension_kwargs),
    ]


setup(
    name = "kg_vect",
    ext_modules = cythonize(ext_modules, annotate=True)
)
