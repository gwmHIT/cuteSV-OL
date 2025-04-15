# coding=utf-8

from setuptools import setup, find_packages
from Cython.Build import cythonize

with open('README.md') as f:
    readme = f.read()

# 使用 Cython 构建所有 .pyx 文件
cython_extensions = cythonize("src/cuteSV/*.pyx")

setup(
    name = "cuteSV-OL",
    version = "1.0.0",
    description = "cuteSV-OL: a real-time structural variation detection framework for nanopore sequencing devices",
    author = "Guo Weimin",
    author_email = "tjiang@hit.edu.cn",
    url = "https://github.com/120L022331/cuteSV-OL",
    license = "MIT",
    packages = find_packages("src"),
    package_dir = {"": "src"},
    package_data={
        "online": ["bin/pandepth"],  # 打包C工具
    },
    data_files = [("", ["LICENSE"])],
    entry_points={
        'console_scripts': [
            'cuteSV=cuteSV.cuteSV:main',
            'cuteSV_ONLINE=online.online:main_function',
        ],
    },
    # long_description = LONG_DESCRIPTION,
    long_description = readme,
    long_description_content_type = 'text/markdown',
    zip_safe = False,
    install_requires = ['scipy', 'pysam', 'Biopython', 'Cigar', 'numpy', 'pyvcf3', 'scikit-learn'],
    ext_modules=cython_extensions
)
