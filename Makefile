#Makefile for developer users
.PHONY: install format lint test sec doc site apidoc prof Cbuild

install:
	pip install -r requirements.txt
format:
	black myfempy
	isort myfempy
lint:
	black myfempy --check
	isort myfempy --check
	flake8 --extend-ignore=E501 myfempy
	interrogate -vv myfempy
test:
	pytest -v
sec:
	pip-audit
doc:
	mkdocs serve
site:
	mkdocs build
apidoc:
	pdoc --output-dir docs/api myfempy/setup/fea.py
prof:
	kernprof -l -v .\__shakedown\profile_test_myfempy.py
#	python -m line_profiler static_line.py.lprof
Cbuild:
	python __dev__setup_wrap_cy_pyx.py build_ext --inplace