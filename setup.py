from setuptools import setup
import re


VERSIONFILE="xfiles/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    version = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))


setup(
    name='xfiles',
    version=version,
    author='',
    author_email='',
    url='',
    description="",
    packages=['xfiles'],
    package_dir={'xfiles': 'xfiles'},
    install_requires=[],
    python_requires='>=3.6',
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    zip_safe=False,
    include_package_data=True,
    classifiers=[],
    project_urls={}
)

