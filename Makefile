#Makefile to dev user
.PHONY: install format lint test sec doc

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