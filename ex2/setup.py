from setuptools import setup, Extension

setup(name='mykmeanssp',
    version='1.0',
    description='k-means C implementation',
    ext_modules=[Extension("mykmeanssp",sources=["kmeans.c"])]
)