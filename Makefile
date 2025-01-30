#Makefile to dev user
.PHONY: install format lint test sec doc perf comp

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
perf:
	kernprof -l -v .\__shakedown\myfempy_perf_test_profile.py
#python -m line_profiler static_line.py.lprof
comp:
	python __dev__setup_wrap_cy_pyx.py build_ext --inplace