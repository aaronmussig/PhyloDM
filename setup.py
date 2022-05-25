import sys

# Restrict to Python 3.7+
if sys.version_info < (3, 7):
    sys.exit('Only Python 3.7+ is supported')

# Check setuptools is installed
try:
    from setuptools import setup
except ImportError:
    sys.exit('Please install setuptools before installing this package.')

# Check setuptools_rust is installed
try:
    from setuptools_rust import Binding, RustExtension
except ImportError:
    sys.exit('Please install setuptools-rust before installing this package.')


# Read the long description from the README file
def readme():
    with open('README.md') as f:
        return f.read()


# Read the version from the Cargo.toml file
def get_package_version():
    with open('Cargo.toml') as f:
        for line in f:
            if line.startswith('version ='):
                return line.split('=')[1].strip().strip('"').strip("'")


# Create the package
setup(name='phylodm',
      version=get_package_version(),
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
          'Development Status :: 4 - Beta',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Natural Language :: English',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 3.7',
          'Programming Language :: Python :: 3.8',
          'Programming Language :: Python :: 3.9',
          'Programming Language :: Python :: 3.10',
          'Programming Language :: Rust',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
      ],
      keywords='phylogenetic distance matrix symmetric',
      packages=['phylodm'],
      install_requires=['numpy'],
      setup_requires=['setuptools-rust', 'setuptools', 'wheel'],
      python_requires='>=3.7',
      data_files=[("", ["LICENSE"])],
      rust_extensions=[RustExtension("phylodm.pdm", binding=Binding.PyO3, features=['python'])],
      zip_safe=False,
      )
