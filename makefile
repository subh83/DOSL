# Discrete Optimal search Library (DOSL)
# Makefile for compiling examples

# Change these if needed or pass as input parameter
PREFIX = /usr/local

# ================================

# INCLUDE FOLDER
# --------------

DOSL_FOLDER = $(PREFIX)/include

# --------------------------------------------
# Install

.PHONY: install
install:
	cp -r dosl $(DOSL_FOLDER)
	find $(DOSL_FOLDER)/dosl -type d -exec chmod 755 {} \;
	find $(DOSL_FOLDER)/dosl -type f -exec chmod 644 {} \;

.PHONY: uninstall
uninstall:
	rm -rf $(DOSL_FOLDER)/dosl

# --------------------------------------------
# Examples

.PHONY: examples
examples:
	cd examples-dosl && make all

.PHONY: examples-test
examples-test:
	cd examples-dosl && make test

.PHONY: examples-simple
examples-simple:
	cd examples-dosl && make simple

.PHONY: examples-advanced
examples-advanced:
	cd examples-dosl && make advanced

.PHONY: examples-run
examples-run:
	cd examples-dosl && make clean

.PHONY: examples-clean
examples-clean:
	cd examples-dosl && make clean

.PHONY: run-examples
run-examples:
	cd examples-dosl && make run

