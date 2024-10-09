Welcome to the NMRInversions.jl package!

Hopefully this documentation will provide all the information you need to get started.

!!! info "Installing julia"
    You can download the julia programming language by following [this link](https://julialang.org/downloads/).
    Afterwards, you can either use the julia console (aka REPL) 
    by finding the julia .exe file in your computer, or you can use VScode,
    which provides an environment similar to MATLAB with the use of the 
    [julia extension](https://www.julia-vscode.org/)
    (follow link for installation instructions). 
    Of course, if you already have a preferred development workflow 
    and you know what you're doing, by all means go for it.

!!! info "Installing the package"
    The package can be installed by running the following command on the julia console:
    ```
    using Pkg ; Pkg.add("NMRInversions")
    ```
    This usually takes a while, but it needs to be done only once 
    (unless you swap environment, more on that 
    [here](https://pkgdocs.julialang.org/v1/environments/)).

    Afterwards, you can use the package by running 
    ```
    using NMRInversions
    ```
    in your julia console, every time you start a new session.
    

    Whenever a new version comes up, you can run:
    ```
    using Pkg ; Pkg.update("NMRInversions")
    ```
    to update the package.

!!! info "GLMakie extension"
    The package provides an extension for interactive visualization,
    using the GLMakie package. To gain access to these capabilities,
    you also need to install GLMakie:
    ```
    using Pkg ; Pkg.add("GLMakie")
    ```
    And to use it, you need to run
    ```
    using GLMakie
    ```
    in order to access the plotting functions.
    

For more details on how to use the package, you can start by refering to the [tutorial](tutorial.md) section.

