import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="squeaky",
    version="0.0.1",
    author="Gerry Tonkin-Hill, Neil MacAlistair, Chris Ruis",
    author_email="g.tonkinhill@gmail.com",
    description="A pan-genome analysis pipeline",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/gtonkinhill/squeaky",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
