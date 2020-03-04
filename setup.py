import setuptools
import os

if os.path.lexists("README.md"):
    fh=open("README.md", "r")
    long_description=fh.read()
    fh.close()
else:
    long_description="`MyBook` Book"

setuptools.setup(
    name='pymcel',  
    author="Author Y. Name",
    author_email="author@email",
    description="Short description of package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/user/repo",
    keywords='keyword here',
    license='MIT',
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        ],
    version='0.0',
    packages=setuptools.find_packages(),
    install_requires=[
        'scipy','ipython','matplotlib',
        ],
    include_package_data=True,
 )
