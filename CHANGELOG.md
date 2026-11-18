# CHANGELOG

All notable changes to **myfempy** (*MultiphYsics Finite Elements Module to PYthon*) will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and [git-changelog](https://github.com/pawamoy/git-changelog), adhering to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## \[0.10.5\] - 2026-10-28

### ✨ Features (`feat`)

* **Core / Cython:** Updated `myfempy` core modules and Cython compilation wrappers (`setup_wrap_cy_pyx.py`) for version `0.10.5`.

---

## \[0.10.4\] - 2026-10-25

### 🐛 Bug Fixes (`fix`)

* **Core / Packaging:** General bug fixes across package modules and `pyproject.toml` configuration.

### 📚 Documentation (`docs`)

* **MkDocs / Examples:** Updated `docs/` documentation and `examples/` scripts for `v0.10.4`.

---

## \[0.10.0\] - 2026-09-15

### ✨ Features (`feat`)

* **Python 3.14 & Free-Threaded (`3.14t`):** Added full support for Python 3.14, including automatic parallelization in the Cython code assembler when running in free-threaded CPython builds.

### 🐛 Bug Fixes (`fix`)

* **Preview Module:** Disabled model preview module under Python `3.14t` to prevent upstream compatibility crashes with `vtk`.

---

[Unreleased]: https://github.com/avgcampos/myfempy/compare/v0.10.5...HEAD
[0.10.5]: https://github.com/avgcampos/myfempy/compare/v0.10.4...v0.10.5
[0.10.4]: https://github.com/avgcampos/myfempy/compare/v0.10.0...v0.10.4
[0.10.0]: https://github.com/avgcampos/myfempy/releases/tag/v0.10.0