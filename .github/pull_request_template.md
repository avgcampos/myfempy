## Summary

<!-- Briefly describe what this Pull Request changes in myfempy and why. -->

## Changelog Category (`git-changelog`)

<!-- Mark the appropriate label(s) with [x] to help categorize the release notes -->

- [ ] 💥 `BREAKING CHANGE`: Incompatible API or solver input changes
- [ ] ✨ `feat`: New feature, finite element, physics module, or solver capability
- [ ] 🐛 `fix`: Bug fix, numerical stability fix, or runtime correction
- [ ] ⚡ `perf`: Performance or parallelization improvement (Cython / Python `3.14t`)
- [ ] ♻️ `refactor`: Code restructuring without changing external behavior
- [ ] 📚 `docs`: Documentation updates (MkDocs, docstrings, or `examples/`)
- [ ] 🎨 `style`: Code formatting (`black`), linting, or CLI output adjustments
- [ ] 🧪 `test`: Unit tests or verification benchmarks (`test_myfempy.py`)
- [ ] 🏗️ `build`: Build system, `pyproject.toml`, `setup_wrap_cy_pyx.py`, or dependencies
- [ ] ⚙️ `ci`: GitHub Actions workflows (`.github/workflows/`)
- [ ] 🔧 `chore`: Routine maintenance, metadata, or repository housekeeping
- [ ] ⏪ `revert`: Reverts a previous commit or feature

## Affected Modules

<!-- Optional: specify the affected areas (e.g., assembler, cython, mesh, post, docs) -->

- 

## Checklist

- [ ] Commit messages follow Conventional Commits (`feat:`, `fix:`, `perf:`, etc.)
- [ ] Code is formatted with `black`
- [ ] Verified locally (including Cython compilation / Python `3.14` & `3.14t` compatibility if applicable)
- [ ] Updated `CHANGELOG.md` under `[Unreleased]` (or ready for automated release notes)