JULIA = julia
RUN_JULIA = time $(JULIA) --color=yes

JULIA_PROJECT = @.
export JULIA_PROJECT

OMP_NUM_THREADS = 2
export OMP_NUM_THREADS

# See [[./.travis.yml::GKS_WSTYPE]]
GKS_WSTYPE ?= png
export GKS_WSTYPE

.PHONY: help test prepare

help:
	@cat misc/make-help.md

Manifest.toml: ci/before_script.jl
	$(RUN_JULIA) $<

prepare:
	$(MAKE) Manifest.toml

test: prepare
	$(RUN_JULIA) --check-bounds=yes test/runtests.jl

include misc/docs.mk
