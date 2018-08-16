JULIA ?= julia
RUN_JULIA = time $(JULIA) --color=yes

# Make variable JULIA is cached in config.mk (which is created by make
# automatically or manually configured by "make configure"):
-include config.mk

JULIA_PROJECT = @.
export JULIA_PROJECT

OMP_NUM_THREADS = 2
export OMP_NUM_THREADS

# See [[./.travis.yml::GKS_WSTYPE]]
GKS_WSTYPE ?= png
export GKS_WSTYPE

.PHONY: help test parallel-test test_* prepare configure

help:
	@cat misc/make-help.md

config.mk:
	echo "JULIA = $(JULIA)" > $@

configure:
	$(MAKE) JULIA=$(JULIA) config.mk --always-make

Manifest.toml: ci/before_script.jl config.mk
	mkdir -pv tmp/Manifest
	-cp -v Manifest.toml tmp/Manifest/backup.$$(date +%Y-%m-%d-%H%M%S).toml
	$(RUN_JULIA) $<

prepare:
	$(MAKE) Manifest.toml

test: prepare
	$(RUN_JULIA) --check-bounds=yes test/runtests.jl

TEST_NAMES = $(patsubst test/test_%.jl, test_%, $(wildcard test/test_*.jl))

parallel-test:
	mkdir -pv tmp/test
	-rm tmp/test/test_*.log
	$(MAKE) $(TEST_NAMES)
	@echo -e '\033[0;32m Finished! \033[0m'

$(TEST_NAMES): test_%: test/test_%.jl prepare
	mkdir -pv tmp/test
	$(RUN_JULIA) --check-bounds=yes $< 2>&1 | tee tmp/test/$@.log

include misc/docs.mk
