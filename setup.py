from setuptools import setup
from setuptools_rust import Binding, RustExtension

setup(
    name="phylotree",
    version="1.0.0",
    rust_extensions=[
        RustExtension("phylotree.pytree", binding=Binding.PyO3, features=["python"])
    ],
    packages=["phylotree"],
    setup_requires=["setuptools-rust", "setuptools", "wheel"],
    zip_safe=False,
)
