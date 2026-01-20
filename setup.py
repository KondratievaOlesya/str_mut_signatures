"""Setup script for str_mut_signatures."""

import re
from pathlib import Path

from setuptools import find_packages, setup

with open('README.rst') as readme_file:
    readme = readme_file.read()

# Read version from __init__.py
init_file = Path(__file__).parent / 'src' / 'str_mut_signatures' / '__init__.py'
with open(init_file) as f:
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", f.read(), re.M)
    if version_match:
        version = version_match.group(1)
    else:
        raise RuntimeError("Unable to find version string.")

requirements = [
    'pandas>=2.0.0',
    'numpy>=1.24',
    'scikit-learn>=1.3',
    'matplotlib>=3.7',
    'trtools'
]

test_requirements = [
    'pytest>=8.0.0',
    'pytest-cov>=4.1.0',
    'pytest-regressions'
]

setup(
    author="Olesia Kondrateva",
    author_email='xkdnoa@gmail.com',
    python_requires='>=3.8',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13',
    ],
    description="Extract STR signatures from annotated VCF",
    entry_points={
        'console_scripts': [
            'str_mut_signatures=str_mut_signatures.cli:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme,
    include_package_data=True,
    keywords='str_mut_signatures',
    name='str_mut_signatures',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/acg-team/str_mut_signatures',
    version=version,
    zip_safe=False,
)
