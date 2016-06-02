import os
import sys
from setuptools import setup
from setuptools.command.test import test as TestCommand


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


def read_requires(fname):
    with open(fname, 'r') as fh:
        requires = [l.strip() for l in fh.readlines()]

    return requires


def read_version():
    for line in open(os.path.join('antismash', '__init__.py'), 'r'):
        if line.startswith('__version__'):
            return line.split('=')[-1].strip().strip("'")


class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        import pytest
        errcode = pytest.main(self.test_args)
        sys.exit(errcode)


setup(
    name='antismash-skeleton',
    version=read_version(),
    author='Kai Blin',
    author_email='kblin@biosustain.dtu.dk',
    description='Skeleton implementation of antiSMASH to test new plugins',
    long_description=read('README.md'),
    install_requires=read_requires('requirements.txt'),
    cmdclass={'test': PyTest},
    entry_points={
        'console_scripts': [
            'skeleton-antismash=antismash.__main__:main'
        ],
    },
    packages=['antismash'],
    url='https://bitbucket.org/kblin/skeleton',
    license='GNU AGPL',
    classifiers=[
        'Programming Language :: Python',
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU Affero General Public License',
        'Operating System :: OS Independent',
    ],
)
