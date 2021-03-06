__author__ = 'junz'

import sys
import io
import os
import codecs
import re
from setuptools import setup, find_packages

here = os.path.abspath(os.path.dirname(__file__))
os.chdir(here)

package_root = 'v1dd_physiology'

# get install requirements
with open('requirements.txt') as req_f:
    install_reqs = req_f.read().splitlines()
install_reqs = [ir for ir in install_reqs if '#' not in ir]
install_reqs = install_reqs[::-1]
print('\ninstall requirements:')
print('\n'.join(install_reqs))
print('')

# get long_description
def read(*filenames, **kwargs):
    encoding = kwargs.get('encoding', 'utf-8')
    sep = kwargs.get('sep', '\n')
    buf = []
    for filename in filenames:
        with io.open(filename, encoding=encoding) as f:
            buf.append(f.read())
    return sep.join(buf)

long_description = read('README.md')

# find version
def find_version(f_path):
    version_file = codecs.open(f_path, 'r').read()
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")
version = find_version(os.path.join(here, 'v1dd_physiology', '__init__.py'))

setup(
    name='v1dd_physiology',
    version=version,
    python_requires='>3.8',
    url='https://github.com/AllenInstitute/v1dd_physiology.git',
    author='Allen Institute',
    install_requires=install_reqs,
    author_email='junz@alleninstitute.org',
    description='SDK for volumetric multi-photon data of V1 DeepDive project',
    long_description=long_description,
    packages=find_packages(),
    include_package_data=True,
    package_data={'':['*.md', '*.txt', '*.hdf5', '.json']},
    platforms='any',
    classifiers=[
        'Programming Language :: Python',
        'Development Status :: 4 - Beta',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        ],
)