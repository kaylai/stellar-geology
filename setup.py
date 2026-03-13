import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="stellar-geology",
    version="0.0.1",
    author="Kayla Iacovino",
    author_email="kaylaiacovino@gmail.com",
    description=("A python library for calculating planet geochemistry from stellar compositions."),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/kaylai/stellar-geology",
    packages=setuptools.find_packages(),
    install_requires=[],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.10',
)
