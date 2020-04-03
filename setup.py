from setuptools import setup
import os
import re


def read_version():
    path = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'phylodm/__init__.py')
    with open(path, 'r') as fh:
        return re.search(r'__version__\s?=\s?[\'"](.+)[\'"]', fh.read()).group(1)


def readme():
    with open('README.md') as f:
        return f.read()


setup(name='phylodm',
      version=read_version(),
      description='Efficient calculation of phylogenetic distance matrices.',
      long_description=readme(),
      long_description_content_type='text/markdown',
      author='Aaron Mussig',
      author_email='aaronmussig@gmail.com',
      url='https://github.com/aaronmussig/PhyloDM',
      license='GPL3',
      project_urls={
          "Bug Tracker": "https://github.com/aaronmussig/PhyloDM/issues",
          "Documentation": "https://github.com/aaronmussig/PhyloDM",
          "Source Code": "https://github.com/aaronmussig/PhyloDM",
      },
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Natural Language :: English',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
          'Programming Language :: Python :: 3.8',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
      ],
      keywords='phylogenetic distance matrix symmetric',
      packages=['phylodm'],
      entry_points={
          'console_scripts': [
              'phylodm = phylodm.__main__:main'
          ]
      },
      install_requires=['numpy', 'dendropy', 'h5py', 'tqdm'],
      python_requires='>=3.6',
      data_files=[("", ["LICENSE"])]
      )
