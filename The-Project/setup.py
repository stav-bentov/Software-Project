from setuptools import setup, Extension

module = Extension(
            "my_spkmeans",
            sources=["spkmeans.c", "kmeans.c", "spkmeansmodule.c"])

setup(name='my_spkmeans',
    version='1.0',
    description='spkmeans C implementation',
    ext_modules=[module]
)