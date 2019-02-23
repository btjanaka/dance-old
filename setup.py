from setuptools import setup


def readme():
    with open("README.md") as f:
        return f.read()


setup(
    name="dance",
    version="0.0.0",  # Matches __init__.py
    author="Bryon Tjanaka",
    author_email="bryon.tjanaka@gmail.com",
    description=("Generate a database of QM calculations for select trivalent "
                 "nitrogen molecules"),
    long_description=readme(),
    long_description_content_type="text/markdown",
    install_requires=[
        "matplotlib==3.0.2",
        "numpy==1.16.1",
        "openeye-toolkits",
    ],
    extras_require={
        "dev": ["pytest", "yapf"],
    },
    license="MIT",
    keywords="chemistry molecules forcefield",
    packages=["dance"],
    scripts=["bin/dance"],
)
