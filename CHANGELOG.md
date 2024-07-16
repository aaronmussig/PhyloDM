# [3.1.0](https://github.com/aaronmussig/PhyloDM/compare/v3.0.0...v3.1.0) (2024-07-16)


### Features

* **edges:** update_edge_lengths and update_all_edge_lengths added (Thanks @FinnOD!) ([bed49be](https://github.com/aaronmussig/PhyloDM/commit/bed49be78bd587cd683fc1ba461c9f922ffd8049))
* **maturin:** Add maturin for building. ([7ab9a46](https://github.com/aaronmussig/PhyloDM/commit/7ab9a4670aa25f4769fc05d136401947586f2b35))
* **maturin:** Add maturin for building. ([010da14](https://github.com/aaronmussig/PhyloDM/commit/010da1401521ebd9dbfc93992cca2d56f0efc64c))
* **maturin:** Add maturin for building. ([64c955b](https://github.com/aaronmussig/PhyloDM/commit/64c955b8d480469e14f904b40baf3bf8556f36c1))
* **rust:** Make all methods public for other applications. ([aaea7e5](https://github.com/aaronmussig/PhyloDM/commit/aaea7e576d2b1dd58bec9069d0e9cbfa18bd953c))

# [3.0.0](https://github.com/aaronmussig/PhyloDM/compare/v2.2.1...v3.0.0) (2023-06-17)


### Features

* **error handling:** Rust Result<> on fallible, Python extended newick support ([4d0520e](https://github.com/aaronmussig/PhyloDM/commit/4d0520e480380adf3f75364159784feec4cf27be))


### BREAKING CHANGES

* **error handling:** for Python.

Rust:
- Breaking API changes for all
fallible methods, Result must be handled.
* **error handling:** Rust API contains Result outputs for all fallible
methods.

## [2.2.1](https://github.com/aaronmussig/PhyloDM/compare/v2.2.0...v2.2.1) (2023-01-09)


### Performance Improvements

* **rust:** Replace leaf node index lookup HashMap with vector to improve performance. ([b857adb](https://github.com/aaronmussig/PhyloDM/commit/b857adbd67538bb58ec75896b0190bfe586ccbdf))
* **rust:** Update iteration methods and reduce overheads when initialising new vectors/maps. ([6fb461a](https://github.com/aaronmussig/PhyloDM/commit/6fb461abdf2df4d66d886bdb27c143aa38ea67f8))

# [2.2.0](https://github.com/aaronmussig/PhyloDM/compare/v2.1.3...v2.2.0) (2023-01-03)


### Features

* **python:** Add distance method to allow for distance calculation without constructing a distance matrix. ([7f86835](https://github.com/aaronmussig/PhyloDM/commit/7f868354dff4f08ecc992517227b072d82b1eff0))

## [2.1.3](https://github.com/aaronmussig/PhyloDM/compare/v2.1.2...v2.1.3) (2022-10-04)


### Bug Fixes

* **bump:** Version bump to rebuild wheels for py37, py38. ([a6a5ade](https://github.com/aaronmussig/PhyloDM/commit/a6a5ade17ac476286b909d1d0c083db32a4891dc))

## [2.1.2](https://github.com/aaronmussig/PhyloDM/compare/v2.1.1...v2.1.2) (2022-07-28)


### Bug Fixes

* **pdm:** Make distance a non-mutable call. ([ebef8e7](https://github.com/aaronmussig/PhyloDM/commit/ebef8e7979666f9d921ffb9641f5f929b4be3da6))

## [2.1.1](https://github.com/aaronmussig/PhyloDM/compare/v2.1.0...v2.1.1) (2022-07-27)


### Bug Fixes

* **python:** Order taxa before calculating distance matrix for reproducibility. ([8dc8d29](https://github.com/aaronmussig/PhyloDM/commit/8dc8d2991fb4998170e609d4e37a82d46459dcef))

# [2.1.0](https://github.com/aaronmussig/PhyloDM/compare/v2.0.5...v2.1.0) (2022-07-27)


### Bug Fixes

* **dm:** Fixed calling DM twice panicking. ([0604030](https://github.com/aaronmussig/PhyloDM/commit/0604030517e54e5e9fb0d580d4d13e91068e8b12))


### Features

* **rust:** Major refactor of API to be more rust-friendly. ([29fe848](https://github.com/aaronmussig/PhyloDM/commit/29fe848ff9fe08393ea08359fe13336ce9e8af86))

## [2.0.5](https://github.com/aaronmussig/PhyloDM/compare/v2.0.4...v2.0.5) (2022-06-01)


### Bug Fixes

* **setup.py:** Fix unicode for setuptools. ([7912966](https://github.com/aaronmussig/PhyloDM/commit/7912966c75c665938daad6d93c2168a75e793138))

## [2.0.4](https://github.com/aaronmussig/PhyloDM/compare/v2.0.3...v2.0.4) (2022-06-01)


### Bug Fixes

* **rust:** Updated to publish to cargo. ([ef32146](https://github.com/aaronmussig/PhyloDM/commit/ef32146be61c94c23bcabbf5f03f6e9794f60b77))

## [2.0.3](https://github.com/aaronmussig/PhyloDM/compare/v2.0.2...v2.0.3) (2022-05-31)


### Bug Fixes

* **ci:** Updated CI for semantic release. ([c5b15dc](https://github.com/aaronmussig/PhyloDM/commit/c5b15dcd8ff6f48c4890079203e24de33d3c2ec0))

## [2.0.2](https://github.com/aaronmussig/PhyloDM/compare/v2.0.1...v2.0.2) (2022-05-31)


### Bug Fixes

* **ci:** Updated CI and fixed build not installing. ([fc1d645](https://github.com/aaronmussig/PhyloDM/commit/fc1d6455d165143b0d2787b9f129e6aefc2221c6))
* **packaging:** Fixed an issue where packaging did not work. ([4ebd2f1](https://github.com/aaronmussig/PhyloDM/commit/4ebd2f1cccc714827d29c33c1ab30c54eaae52d2))
