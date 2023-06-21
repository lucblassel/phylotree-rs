from setuptools import setup
from setuptools_rust import Binding, RustExtension
import tomllib as toml


def version():
    with open("Cargo.toml", "rb") as file:
        cargo = toml.load(file)
    return cargo["package"]["version"]


def readme():
    with open("README.md") as file:
        return file.read()


setup(
    name="phylotree",
    version=version(),
    description="A Rust backed package to deal with phylogenetic trees",
    long_description=readme(),
    long_description_content_type="text/markdown",
    author="Luc Blassel",
    license="GPL3",
    project_urls={
        "Bug Tracker": "https://github.com/lucblassel/phylotree-rs/issues",
        "Documentation": "https://github.com/lucblassel/phylotree-rs",
        "Source Code": "https://github.com/lucblassel/phylotree-rs",
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Rust",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.7",
    keywords="phylogenetic tree distance matrix",
    data_files=[("", ["LICENSE"])],
    rust_extensions=[
        RustExtension("phylotree.pytree", binding=Binding.PyO3, features=["python"])
    ],
    packages=["phylotree"],
    setup_requires=["setuptools-rust", "setuptools", "wheel"],
    zip_safe=False,
)
