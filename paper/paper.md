---
title: 'NMRInversions.jl, a julia package for time-domain Nuclear Magnetic Resonance.'
date: xx/xx/xxxx 
tags:
  - Magnetic Resonance
  - NMR relaxation
  - NMR diffusion
  - julia
  - Inverse problems
  - Numerical inversion
authors:
  - name: Aristarchos Mavridis
    affiliation: "1"
    orcid: 0000-0002-6619-2303
  - name: Carmine D'Agostino
    affiliation: "1,2"
    orcid: 0000-0003-3391-8320

affiliations:
  - index: 1
    name: Department of Chemical Engineering, The University of Manchester, Oxford Road, Manchester, UK
  - index: 2
    name: Dipartimento di Ingegneria Civile, Chimica, Ambientale e dei Materiali (DICAM), Alma Mater Studiorum – Università di Bologna, Via Terracini, 28, 40131 Bologna, Italy

bibliography: paper.bib
---

https://joss.readthedocs.io/en/latest/example_paper.html

# Summary

This package provides a library of functions as a user-friendly interface, 
to perform time-domain NMR data processing and visualizations.
It is aimed towards people who are not necessarily familiar with the julia programming language.

Functionality includes importing data from various NMR instrument formats,
performing phase-correction on raw data, fitting exponential decays using numerical inversions, and visualizing these results.

# Statement of need

NMR relaxation and diffusion methods are popular in various fields of science 
and engineering, with applications ranging from studying the properties of 
porous rocks, to catalyst supports, biological tissues and many more.

However, a lot of the users of such experimental techniques are not fully familliar with
the technical aspects of the numerical inversions involved in the data processing, 
and sometimes find it challenging to write their own algorithms.

This package can cover the needs of users who would like a simple interface for everyday data processing, 
while retaining the ability to dive into the details and have absolute control over every step of the process.

The julia programming language [@bezanson2017julia] is an excellent option for such a package, since it 
provides user-friendly, high level syntax which is similar to MATLAB, whitout any 
compromise when it comes to computational performance and scalability. It's a thriving ecosystem for scientific
computing applications, with several optimization libaries able to tackle the problems arising in NMR applications. 

This can enable easy integration with the latest advances from the literature, 
without the need for writing solvers in low-level languages for performance purposes, 
in line with julia's philosophy of solving the "two language problem".

Additionally, julia offers multiple dispatch capabilities, which allows us to use the same function names for different operations,
depending on the function inputs.
This enables our package to be more user-friendly, as it really reduces the ammount of function names to be memorized.

Instrument manufactures often provide algorithms similar to the ones provided in this library, however, 
they are often part of GUI interfaces, which don't offer costumization opportunities, 
and lack the flexibility of writing scripts to efficiently process large sets of data. 

Furthermore, general efforts should be made in improving the reproducibility of scientific research, 
and having accessible, open-source algorithms is of paramount importance when working towards such goals.

# Theory

Numerical inversion methods have been used extensively in the NMR literature for the last fifty years.
Recent advances in computer performance have made these methods more accessible than ever, 
with the average laptop nowadays being easily capable of tackling such computations with ease.

The methods used in this study are mainly based on Michell et. al. 2012, 
and the mathematical notation we use in the source code mostly follows the aforementioned paper. 

Relaxation and diffusion NMR experiments in fluids yield exponentially-decaying signals, according to the BPP theory.
Characterizing the rate of these decays can provide valuable information on the physical and 
chemical characteristics of the sample.

Usually, there are various relaxation times and/or diffusion coefficients within most samples of interest, 
such as different pore environments, chemical sites and magnetic field inhomogeneities, 
and therefore multi-exponential forms are commonly used to characterize the data.

For example, the magnetization signal ($M$) from a CPMG experiment can be described using the equation:

$$ M = \sum_{n=1}^{N} \alpha_n exp({t/T_2_n}) $$

For small numbers of N, (1 to 3) finding the least squares solution for  the $\alpha$ and $T_2$ 
values is a simple nonlinear regression problem.

However, we are often interested in distributions of many relaxation times, 
e.g., corresponding to pore size distributions in rock samples or other materials.

The expression above can also written using an almost equivalent, matrix multiplication form:

$$ K*f = g ,$$

where $g$ is a vector of experimental results, $f$ contains the $/alpha$ values, and $K$ is the kernel matix, 
whose columns are exponential decays of several relaxation times within a predetermined range (often logarithmically spaced).

One might be tempted to use the inverse of $K$ and solve the problem as $ f = K^{-1}g $ , but this would 
produce very unreliable results, due to the high condition number of the kernel making the solution extremely sensitive to the noise in $g$.

A common method to obtain more reliable solutions to use Tikhonov regularization, and solve the following minimization problem instead.
$$ \min_{f \geq 0} ||Kf-g||_2^2 + \alpha||f||_2^2  $$


# Example use

# Acknowledgements
The authors would like to acknowledge funding from bp-ICAM and the EPSRC grant no. EP/V519613/1. 
Furthermore, A. Mavridis would like to thank the julia community for numerous helpful discussions around the topic.

# References
