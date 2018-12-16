init:
	pip install -r requirements.txt

setup:
	setup.py

test:
	py.test tests

.PHONY: init test