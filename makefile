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

.PHONY: examples
examples:
	cd examples-dosl && make all

.PHONY: uninstall
uninstall:
	rm -rf $(DOSL_FOLDER)/dosl
