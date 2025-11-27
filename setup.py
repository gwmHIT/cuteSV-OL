# coding=utf-8

from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize

with open('README.md') as f:
    readme = f.read()

# 使用 Cython 构建所有 .pyx 文件
cython_extensions = cythonize([
    Extension(
        "cuteSV.genotype_improve",                  # 模块名：包名.文件名
        ["src/cuteSV/genotype_improve.pyx"],        # 对应的 pyx 源文件
        extra_compile_args=["-std=c99"],            # 关键：开启 C99
    ),
    Extension(
        "cuteSV.prove_func",
        ["src/cuteSV/prove_func.pyx"],
        extra_compile_args=["-std=c99"],
    ),
])

setup(
    name = "cuteSV-OL",
    version = "1.0.0",
    description = "cuteSV-OL: a real-time structural variation detection framework for nanopore sequencing devices",
    author = "Guo Weimin",
    author_email = "tjiang@hit.edu.cn",
    url = "https://github.com/gwmHIT/cuteSV-OL",
    license = "MIT",
    packages = find_packages("src"),
    package_dir = {"": "src"},
    package_data={
        "online": ["bin/*", "test/*"],
    },
    data_files = [("", ["LICENSE"])],
    entry_points={
        'console_scripts': [
            'cuteSV_RT=cuteSV.cuteSV:main',
            'cuteSV_ONLINE=online.online:main_function',
        ],
    },
    # long_description = LONG_DESCRIPTION,
    long_description = readme,
    long_description_content_type = 'text/markdown',
    zip_safe = False,
    ext_modules=cython_extensions
)
