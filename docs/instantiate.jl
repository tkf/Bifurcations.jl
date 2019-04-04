#!/bin/bash
# -*- mode: julia -*-
#=
JULIA="${JULIA:-julia --color=yes --startup-file=no}"
export JULIA_PROJECT="$(dirname ${BASH_SOURCE[0]})"
exec ${JULIA} "${BASH_SOURCE[0]}" "$@"
=#

using Pkg
Pkg.instantiate()
Pkg.develop(PackageSpec(path=dirname(@__DIR__)))
