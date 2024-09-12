Here are some examples on how to use the package.

# 1D inversion
Suppose we would like to perform an inversion for a CPMG experiment.
We need some data to work with, so let's create some.

```julia
using NMRInversions

# Define the experiment's time acquisition range (in seconds)
x = collect(range(0.0001,2,500)) 

# Range of relaxation times
X = exp10.(range(-5, 1, 128)) 

# Create a kernel matrix
K = create_kernel(CPMG, x, X) 

# Create a distribution of relaxation times
f_custom = [0.5exp.(-(x)^2 / 3) + exp.(-(x - 1.3)^2 / 0.5) for x in range(-5, 5, length(X))]

# Use the "forward problem" to create the NMR signal from the T2 distribution
g = K * f_custom

# Add some noise to the data
y = g + 0.005 * maximum(g) .* randn(length(x))
```

Now we have the data, the inversion can be performed using a single line of code!
```julia
results = invert(CPMG, x, y, alpha=gcv)
```

And the results can easily be visualised using the GLMakie extension of the package.
```julia
using GLMakie
plot(results)
```

The `results` variable is an `inv_out_1D` structure, which contains all the relevant information produced by the `inversion` function.
To access that information, we can look at the fields of the structure using the dot notation.
The field names can be shown by using the REPL help mode (typing ? at the julia> prompt), 
and typing the variable's name (in this case ?results). 
Alternatively, running `@doc results` will also give you the same answers.
For example, the vector of values for the relaxation time distribution can be accessed as `results.f`
and the regularization parameter can be accessed as `results.alpha`.

In the example above, the arguments we fed into the `invert()` function were the pulse sequence, 
the time values, and the NMR signal corresponding to each time value. 
Importing data using function such as `import_spinsolve` make things simpler, 
since they return an `input1D` structure which can be used directly in the `invert` function as:
```julia
invert(import_spinsolve())
```

# 2D inversion
Let's try with a T1-T2 experiment, using the IR-CPMG sequence.
Again, we start by creating some data.
