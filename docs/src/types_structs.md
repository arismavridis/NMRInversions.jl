
# Types and Structures

In order to take advantage of Julia's multiple dispatch, 
we need to define some structures which can be used as input to our functions.

Structures are object which may contain multiple fields.
The fields can be accessed using the dot syntax.
For example, if we have a structure named `foo`, with two fields, `a`, and `b`, 
we can access the value of `a` using `foo.a` and the value of `b` using `foo.b`.


## Pulse sequences
```@docs
IR
SR
CPMG
PFG
IRCPMG
```

## Inversion solvers
```@docs
brd
pdhgm
optim_nnls
```

## Finding optimal alpha
```@docs
gcv
lcurve
```

## Inversion and expfit inputs
```@docs
input1D
input2D
```

##  Inversion and expfit outputs
```@docs
inv_out_1D
inv_out_2D
expfit_struct
```


## Kernel structure

```@docs
svd_kernel_struct
```
