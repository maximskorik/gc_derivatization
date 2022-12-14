# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.1]

### Added
* documented functions with docstrings to make use of Python's built-in `help()` function [#23](https://github.com/RECETOX/gc-meox-tms/pull/23)
* added `Publish to PyPi` GitHub Actions workflow [#24](https://github.com/RECETOX/gc-meox-tms/pull/24)

## [1.0.0]

### Added
* added **Anaconda build**, **Python Package with pip** (inc. **SonarCloud**), and **Python Package with Conda** GH Actions.
* added test coverage for main and IO functionality
* added conda dev environment and conda meta.yaml recipe
* added CHANGELOG.md

### Changed
* changed package structure
* divided main functionality, IO handling, and CLI into designated modules

### Fixed
* fixed not working examples in IPython notebook
