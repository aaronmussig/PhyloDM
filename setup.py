import os
import platform
import re
from distutils.command.build import build as build_orig

from setuptools import setup, Extension


def read_version():
    path = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'phylodm/__init__.py')
    with open(path) as fh:
        return re.search(r'__version__\s?=\s?[\'"](.+)[\'"]', fh.read()).group(1)


def readme():
    with open('README.md') as f:
        return f.read()


class build(build_orig):

    def finalize_options(self):
        super().finalize_options()
        __builtins__.__NUMPY_SETUP__ = False
        import numpy
        for extension in self.distribution.ext_modules:
            extension.include_dirs.append(numpy.get_include())
        from Cython.Build import cythonize
        self.distribution.ext_modules = cythonize(self.distribution.ext_modules,
                                                  language_level=3)


compile_extra_args = ['-O3']
link_extra_args = list()
# if platform.system() == "Windows":
#     compile_extra_args = ["/std:c++latest", "/EHsc"]
if platform.system() == "Darwin":
    compile_extra_args.extend(['-std=c++11', "-mmacosx-version-min=10.9"])
    link_extra_args.extend(["-stdlib=libc++", "-mmacosx-version-min=10.9"])

ext_modules = [Extension('phylodm.pdm_c', ['phylodm/pdm_c.pyx'],
                         language='c++',
                         extra_compile_args=compile_extra_args,
                         extra_link_args=link_extra_args
                         )]

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
      package_data={'phylodm': ['pdm_c.pyx']},
      entry_points={
          'console_scripts': [
              'phylodm = phylodm.__main__:main'
          ]
      },
      install_requires=['numpy', 'dendropy', 'h5py', 'tqdm'],
      setup_requires=['cython', 'numpy'],
      python_requires='>=3.6',
      data_files=[("", ["LICENSE"])],
      ext_modules=ext_modules,
      cmdclass={'build': build}
      )
