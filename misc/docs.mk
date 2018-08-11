DOCS_PORT = 35730
LIVERELOAD = livereload --port $(DOCS_PORT) --wait 0.5

.PHONY: docs-* only-serve

docs-build: prepare
	$(RUN_JULIA) docs/build.jl
	cd docs && mkdocs build

docs-cont-build:
	rg --files | entr $(MAKE) docs-build

only-serve:
	cd docs/site && $(LIVERELOAD)
#	cd docs && mkdocs serve --dev-addr=localhost:$(DOCS_PORT)
# Somehow, "mkdocs serve" version didn't work well; it didn't update
# the page on-the-fly.  Use livereload to workaround this.

docs-serve:
	$(MAKE) docs-cont-build &
	$(MAKE) only-serve
