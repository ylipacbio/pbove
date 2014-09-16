SHELL = /bin/bash -e

all: build install

build:
	python setup.py build --executable="/usr/bin/env python"

bdist:
	python setup.py build --executable="/usr/bin/env python"
	python setup.py bdist --formats=egg

install:
	python setup.py install

develop:
	python setup.py develop

test:
	# Unit tests
	find tests/unit -name "*.py" | xargs nosetests
	# End-to-end tests
	@echo pbove cram tests require blasr installed.
	find tests/cram -name "*.t" | xargs cram 

doc:
	sphinx-apidoc -T -f -o doc src/ && cd doc && make html

docs: doc

doc-clean:
	rm -f doc/*.html

clean: doc-clean
	rm -rf dist/ build/ *.egg-info
	rm -rf doc/_build
	find . -name "*.pyc" | xargs rm -f
	rm -rf dist/
	rm -f nostests.xml

pip-install: 
	@which pip > /dev/null
	@pip freeze|grep 'pbtools.pbove=='>/dev/null \
      && pip uninstall -y pbtools.pbove \
      || echo -n ''
	@pip freeze|grep 'pbove=='>/dev/null \
      && pip uninstall -y pbove \
      || echo -n ''
	@pip install --no-index \
          --install-option="--install-scripts=$(PREFIX)/bin" \
          ./

.PHONY: all build bdist install develop test doc clean
