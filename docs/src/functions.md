This page contains the documentation for various useful functions that can be used.

From the julia command line, you can enter '?', 
followed by the name of any function, struct, 
or object you want to learn more about (try it!).

In julia, function definitions look like this:
`foo(x, y ; z)`
Where `x` and `y` are the necessary arguments for the function to work, 
and everything after the semicolon `;` is an optional keyword argument.
In the example above, we can call the function by typing `foo(1, 2)`,
(in which case the default value for `z` will be used), or by typing `foo(1, 2, z=3)`,
to specify the value for z. Sometimes, if there are many keyword arguments, we write 
the function as foo(x, y ; kwargs...). 
For the necessary arguments before the semicolon, order matters. 
For the keyword arguments after the semicolon, order does not matter, but the name of each argument must be specified.

# Importing data functions

# Inversion functions
The most important function is `invert()`, which is the main function of the package.
It works as follows:

```@docs
invert(::Type{<:pulse_sequence1D}, ::AbstractArray, ::Vector)
invert(::input1D)
```

Due to julia's multiple dispatch, 
it is possible to define a function with the same name
but different arguments, to achieve different results.


Because of that, the inversion function also works for 2D inversions,
if the following arguments are used instead:

```@docs
invert(::Type{<:pulse_sequence2D}, ::AbstractVector, ::AbstractVector, ::AbstractMatrix)
invert(::input2D)
```


# Exponential fit functions

For 1D data, we can use the `expfit` function to perform multiexponential fits.
We can use the function by specifying either the number of exponential components,
or a vector which defines the starting points for the regression.
See below:

```@docs
expfit(::Int, ::Type{<:NMRInversions.pulse_sequence1D}, ::Vector, ::Vector )
expfit(::Int, ::input1D)
```  
  
If you have some rough clues about what results you expect, 
it's best to define some starting points close to these, as shown below.
(especially important if you're using double or tri-exponential fits.)

```@docs
expfit(::Vector{Real}, ::Type{<:NMRInversions.pulse_sequence1D}, ::Vector, ::Vector )
expfit(::Vector{Real}, ::input1D)

```


# Plotting functions
This package offers plotting capabilities, using its GLMakie extension.
Simply use `using GLMakie` to load the package, and these functions become available.

We basically take the `plot()` function offered by GLMakie and extend it to types from the package.

For example, for 1D inversions, we have:

```@docs
plot(::NMRInversions.inv_out_1D)
plot!(::Union{Makie.Figure,Makie.GridPosition}, ::NMRInversions.inv_out_1D )

```

For expfits, we have:

```@docs
plot(::NMRInversions.expfit_struct)
plot!(::Union{Makie.Figure,Makie.GridPosition}, ::NMRInversions.expfit_struct )
```

And for 2D inversions, there's an interactive gui to characterize the inversion resutls:

```@docs
plot(::NMRInversions.inv_out_2D)
plot(::NMRInversions.inv_out_2D,::String)
plot!(::Union{Makie.Figure,Makie.GridPosition}, ::NMRInversions.inv_out_2D )
```

