from setuptools import setup, Extension

setup(name='my_spkmeans',
    version='1.0',
    description='k-means C implementation',
    ext_modules=[Extension("my_spkmeans",sources=["spkmeansmodule.c"])]
)