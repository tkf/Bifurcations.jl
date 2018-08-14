    make test

Run `test/runtests.jl`.

    make parallel-test [-j N]

Run tests in parallel with `N` processes.

    make prepare [--always-make]

Run `Pkg.build()` to install required dependencies.
It also does a bit of manual package installation at the moment.

    make configure JULIA={julia Executable}

Configure the julia executable to be used (by default).  If this
command is not run, `julia` found in the `$PATH` is used.

    make docs-serve

Build document and serve it in background.

    make docs-build

Build document by pre-processing Markdown files using Documenter.jl
and run `mkdocs` to convert them to HTML files.
