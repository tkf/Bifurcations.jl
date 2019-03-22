DOCS_PORT = 35730
LIVERELOAD = livereload --port $(DOCS_PORT) --wait 0.5

.PHONY: docs-* only-serve

docs-build: prepare
	$(RUN_JULIA) docs/build.jl

docs-cont-build:
	rg --files | entr $(MAKE) docs-build

only-serve:
	cd docs/build && $(LIVERELOAD)

docs-serve:
	$(MAKE) docs-cont-build &
	$(MAKE) only-serve
