var documenterSearchIndex = {"docs":
[{"location":"functions/#Intro","page":"Functions","title":"Intro","text":"","category":"section"},{"location":"functions/","page":"Functions","title":"Functions","text":"This page contains the documentation for various useful  functions in the NMRInversions package.","category":"page"},{"location":"functions/","page":"Functions","title":"Functions","text":"tip: Tip\nFrom the julia command line, you can enter '?',  followed by the name of any function, struct,  or object you want to learn more about (try it!). After typing using NMRInversions in the julia console,  this feature will work for all the functions mentioned below.","category":"page"},{"location":"functions/","page":"Functions","title":"Functions","text":"info: Info\nIn julia, function definitions look like this:foo(x, y, z ; a, b)For the example above, foo is the name of the function,  and the contents of the parenteses are the arguments.   Within the parenteses, we got two types of arguments:Positional argumentsEverything that appears before the semicolon ;  (x,y and z in this example) is necessary, and must be given in a specific order.Keyword arguments.Everything that appears after the semicolon ;  (b and a in this example) is optional, and can be given in any order, but its name must be specified.In the example above, we can call the function by typing foo(1, 2, 3), (in which case, x=1, y=2, z=3, and the default value for a and b will be used).  You can also call the function by typing foo(1, 2, 3, a=3), to specify the value for a,  or by typing foo(1, 2, 3, b=3, a = 2), to specify the value for both a and b.Sometimes, if there are many keyword arguments, we write  the function as foo(x, y ; kwargs...).  For the necessary arguments before the semicolon, order matters.  For the keyword arguments after the semicolon, order does not matter,  but the name of each argument must be specified. For more information, please refer to this link.","category":"page"},{"location":"functions/#Importing-data-functions","page":"Functions","title":"Importing data functions","text":"","category":"section"},{"location":"functions/","page":"Functions","title":"Functions","text":"This package offers some functions to import NMR experiment data of various formats. If a format you're working with is not yet supported,  please submit an issue  and we'll work on it.","category":"page"},{"location":"functions/","page":"Functions","title":"Functions","text":"The most basic use case would be using data saved in a csv format,  where there are only two columns,  one for your x-axis (time for relaxation and b-factor for diffusion) and one for your y-axis (signal intensity).","category":"page"},{"location":"functions/","page":"Functions","title":"Functions","text":"import_csv(::Type{<:pulse_sequence1D},::String)","category":"page"},{"location":"functions/#NMRInversions.import_csv-Tuple{Type{<:pulse_sequence1D}, String}","page":"Functions","title":"NMRInversions.import_csv","text":"import_csv(seq, file)\n\nImport data from a CSV file. The function reads the file and returns an input1D structure.\n\nseq is the 1D pulse sequence (e.g. IR, CPMG, PFG)\nfile is the path to the CSV file which contains the data (x, y) in two respective columns.\n\nThe function can be called without the seq argument, and the output will be the x and y vectors  (use it like, x,y =import_csv()). Alternatively, the function can also be called with only the seq argument, in which case a file dialog will open to select the file (use it like, data = import_csv(IR))\n\nPlease note that this function will just import the data as is, without any unit conversions. Ensure that your x-axis is in SI.\n\n\n\n\n\n","category":"method"},{"location":"functions/","page":"Functions","title":"Functions","text":"If you're using a spinsolve instrument, you can use the import_spinsolve function. This one requires two files as an input.  The aqcu.par is automatically exported by SpinsolveExpert,  but you might have to export your data file manually in a csv format.","category":"page"},{"location":"functions/","page":"Functions","title":"Functions","text":"import_spinsolve(::String)","category":"page"},{"location":"functions/#NMRInversions.import_spinsolve-Tuple{String}","page":"Functions","title":"NMRInversions.import_spinsolve","text":"import_spinsolve(files)\n\nImport data from a Spinsolve experiment.  Two paths must be provided as follows (order is not important):\n\nfiles = [.../datafile.csv , .../acqu.par] \n\nCalling this function without an argument by typing import_spinsolve() will open a file dialog to select the files. The function reads the acqu.par.bak file to get the acquisition parameters, and the .dat file to get the data.  The function returns an input2D structure.\n\n\n\n\n\n","category":"method"},{"location":"functions/","page":"Functions","title":"Functions","text":"For geospec instruments, you can export your raw data as a text file. That text file can be read by the import_geospec function.","category":"page"},{"location":"functions/","page":"Functions","title":"Functions","text":"import_geospec(::String)","category":"page"},{"location":"functions/#NMRInversions.import_geospec-Tuple{String}","page":"Functions","title":"NMRInversions.import_geospec","text":"import_geospec(dir)\n\nImport data from a .txt format, as exported by Geospec instruments.\n\nThe function reads the relevant information, performs a phase correction on the data, and returns an input1D or input2D structure.\n\nCalling this function without an argument by typing import_geospec()  will open a file dialog to select the .txt file.\n\n\n\n\n\n","category":"method"},{"location":"functions/#Inversion-functions","page":"Functions","title":"Inversion functions","text":"","category":"section"},{"location":"functions/","page":"Functions","title":"Functions","text":"The most important function is invert(), which is the main function of the package. It works as follows:","category":"page"},{"location":"functions/","page":"Functions","title":"Functions","text":"invert(::Type{<:pulse_sequence1D}, ::AbstractArray, ::Vector)\ninvert(::input1D)","category":"page"},{"location":"functions/#NMRInversions.invert-Tuple{Type{<:pulse_sequence1D}, AbstractArray, Vector}","page":"Functions","title":"NMRInversions.invert","text":"Inversion for 1D pulse sequences:\n\ninvert(seq, x, y ; lims, alpha, solver)\n\nThis function will build a kernel and use it to perform an inversion using the algorithm of your choice. The output is an inv_out_1D structure.\n\nNecessary (positional) arguments:\n\nseq is the 1D pulse sequence (e.g. IR, CPMG, PFG)\nx is the experiment x axis (time or b factor etc.)\ny is the experiment y axis (intensity of the NMR signal)\n\nOptional (keyword) arguments:\n\nlims=(a,b,c) will set the \"limits\" of the output X, \n\nso that it starts from 10^a, ends in 10^b and consists of c  logarithmically spaced values  (default is (-5, 1, 128) for relaxation and (-10, -7, 128) for diffusion).  Alternatiively, a vector of values can be used directly, if more freedom is needed  (e.g. lims=exp10.(range(-5, 1, 128))).\n\nalpha determines the smoothing term. Use a real number for a fixed alpha.  No selection will lead to automatically determining alpha through the defeault method, which is gcv.\nsolver is the algorithm used to do the inversion math. Default is brd.\n\n\n\n\n\n","category":"method"},{"location":"functions/#NMRInversions.invert-Tuple{input1D}","page":"Functions","title":"NMRInversions.invert","text":"invert(data::input1D ; kwargs...)\n\nInstead of the positional arguments seq, x and y, you can use a single input1D structure, which contains the same information.  Especially useful if you're using the output of one  of the import functions (look documentation tutorial section).\n\n\n\n\n\n","category":"method"},{"location":"functions/","page":"Functions","title":"Functions","text":"Due to julia's multiple dispatch,  it is possible to define a function with the same name but different arguments, to achieve different results.","category":"page"},{"location":"functions/","page":"Functions","title":"Functions","text":"Because of that, the inversion function also works for 2D inversions, if the following arguments are used instead:","category":"page"},{"location":"functions/","page":"Functions","title":"Functions","text":"invert(::Type{<:pulse_sequence2D}, ::AbstractVector, ::AbstractVector, ::AbstractMatrix)\ninvert(::input2D)","category":"page"},{"location":"functions/#NMRInversions.invert-Tuple{Type{<:pulse_sequence2D}, AbstractVector, AbstractVector, AbstractMatrix}","page":"Functions","title":"NMRInversions.invert","text":"Inversion for 2D pulse sequences:\n\ninvert(seq, x_direct, x_indirect, X_direct, X_indirect, Data)\n\nThis function will build a kernel and use it to perform an inversion using the algorithm of your choice. The output is an inv_out_2D structure.\n\nNecessary (positional) arguments:\n\nseq is the 2D pulse sequence (e.g. IRCPMG)\nx_direct is the direct dimension acquisition parameter (e.g. the times when you aquire CPMG echoes).\nx_indirect is the indirect dimension acquisition parameter (e.g. all the delay times τ in your IR sequence).\nData is the 2D data matrix of complex data.\n\nOptional (keyword) arguments:\n\nlims1 determines the output \"range\" of the inversion in the direct dimension (e.g. T₂ times in IRCPMG)\nlims2 determines the output \"range\" of the inversion in the indirect dimension (e.g. T₁ times in IRCPMG)\n\nIn both cases above, you can use a tuple specifying the limits of the range, or a vector of values, same as the lims argument in the 1D inversion.\n\nalpha determines the smoothing term. Use a real number for a fixed alpha.  No selection will lead to automatically determining alpha through the defeault method, which is gcv.\nsolver is the algorithm used to do the inversion math. Default is brd.\n\n\n\n\n\n","category":"method"},{"location":"functions/#NMRInversions.invert-Tuple{input2D}","page":"Functions","title":"NMRInversions.invert","text":"invert(data::input2D ; kwargs...)\n\nInstead of the positional arguments seq, x_direct , x_indirect and Data, you can use a single input2D  structure, which contains the same information.  Especially useful if you're using the output of one of the  import functions (look documentation tutorial section).\n\n\n\n\n\n","category":"method"},{"location":"functions/#Exponential-fit-functions","page":"Functions","title":"Exponential fit functions","text":"","category":"section"},{"location":"functions/","page":"Functions","title":"Functions","text":"For 1D data, we can use the expfit function to perform multiexponential fits. We can use the function by specifying either the number of exponential components, or a vector which defines the starting points for the regression. See below:","category":"page"},{"location":"functions/","page":"Functions","title":"Functions","text":"expfit(::Int, ::Type{<:NMRInversions.pulse_sequence1D}, ::Vector, ::Vector )\nexpfit(::Int, ::input1D)","category":"page"},{"location":"functions/#NMRInversions.expfit-Tuple{Int64, Type{<:pulse_sequence1D}, Vector, Vector}","page":"Functions","title":"NMRInversions.expfit","text":"expfit(n, seq, x, y;solver)\n\nFit an exponential function to the data x, y using n exponential terms. \n\nStarting points for the nonlinear regression are automatically chosen. The outut is an expfit_struct structure.\n\nArguments:\n\nn : number of exponential terms.\nseq : pulse sequence.\nx : acquisition x parameter (time for relaxation or b-factor for diffusion).\ny : acquisition y parameter (magnetization).\nsolver : optimization solver (default is BFGS).\n\nNote that initial conditions are automatically determined.  If you get NaN for any of the resulting parameters,  try changing the initial conditions by calling the funtion using the u0 argument instead of n.\n\n\n\n\n\n","category":"method"},{"location":"functions/#NMRInversions.expfit-Tuple{Int64, input1D}","page":"Functions","title":"NMRInversions.expfit","text":"expfit(u0, data::input1D; kwargs...)\n\nSimilar to the invert fucntion, expfit can be called using an input1D structure.\n\n\n\n\n\n","category":"method"},{"location":"functions/","page":"Functions","title":"Functions","text":"If you have some rough clues about what results you expect,  it's best to define some starting points close to these, as shown below. (especially important if you're using double or tri-exponential fits.)","category":"page"},{"location":"functions/","page":"Functions","title":"Functions","text":"expfit(::Vector{Real}, ::Type{<:NMRInversions.pulse_sequence1D}, ::Vector, ::Vector )\nexpfit(::Vector{Real}, ::input1D)\n","category":"page"},{"location":"functions/#NMRInversions.expfit-Tuple{Vector{Real}, Type{<:pulse_sequence1D}, Vector, Vector}","page":"Functions","title":"NMRInversions.expfit","text":"expfit(u0, seq, x, y; solver=BFGS())\n\nFit an exponential function to the data x, y. \n\nThe outut is an expfit_struct structure.\n\nArguments:\n\nu0 : initial parameter guesses.\nseq : pulse sequence.\nx : acquisition x parameter (time for relaxation or b-factor for diffusion).\ny : acquisition y parameter (magnetization).\nsolver : OptimizationOptimJL solver, defeault choice is BFGS().\n\nThe u0 argument is a vector of initial parameter guesses,  and it also determines the number of exponential terms used. It should be of the form [a1, b1, a2, b2, ...],  where a's are the amplitudes and b's are the reciprocals of the decay constants. The length of u0 must be an even number.\n\nThe following examples might help to clarify: \n\nexpfit([a,b] , CPMG, x, y) -> mono-exponential fit with initial guess: a * exp.( (1/b) * x) \n\nexpfit([a,b,c,d] , CPMG, x, y) -> double-exponential fit with initial guess: a * exp.( (1/b) * x) + c * exp.((1/d) * x) \n\nexpfit([a,b,c,d,e,f] , CPMG, x, y) -> triple-exponential fit with initial guess: a * exp.( (1/b) * x) + c * exp.((1/d) * x) + e * exp.((1/f) * x) \n\n(where a,b,c,d,e,f are numbers of your choice) Numbers of parameters beyond tri-exponential can also be used, but it is not recommended.\n\n\n\n\n\n","category":"method"},{"location":"functions/#NMRInversions.expfit-Tuple{Vector{Real}, input1D}","page":"Functions","title":"NMRInversions.expfit","text":"expfit(n, data::input1D; kwargs...)\n\nSimilar to the invert fucntion, expfit can be called using an input1D structure.\n\n\n\n\n\n","category":"method"},{"location":"functions/#Plotting-functions","page":"Functions","title":"Plotting functions","text":"","category":"section"},{"location":"functions/","page":"Functions","title":"Functions","text":"This package offers plotting capabilities, using its GLMakie extension. Simply use using GLMakie to load the package, and these functions become available.","category":"page"},{"location":"functions/","page":"Functions","title":"Functions","text":"We basically take the plot() function offered by GLMakie and extend it to types from the package.","category":"page"},{"location":"functions/","page":"Functions","title":"Functions","text":"For example, for 1D inversions, we have:","category":"page"},{"location":"functions/","page":"Functions","title":"Functions","text":"plot(::NMRInversions.inv_out_1D)\nplot!(::Union{Makie.Figure,Makie.GridPosition}, ::NMRInversions.inv_out_1D )\n","category":"page"},{"location":"functions/#MakieCore.plot-Tuple{inv_out_1D}","page":"Functions","title":"MakieCore.plot","text":"plot(res::NMRInversions.inv_out_1D...; kwargs...)\n\nPlot the results contained in a inv_out_1D structure. This function can take any number of inv_out_1D structures as input. \n\n\n\n\n\n","category":"method"},{"location":"functions/#MakieCore.plot!-Tuple{Union{GridPosition, Figure}, inv_out_1D}","page":"Functions","title":"MakieCore.plot!","text":"plot!(fig, res...; title)\n\nPlot the results contained in a inv_out_1D structure on a figure or a grid position object. This function can take any number of inv_out_1D structures as input.\n\nThe arguments are:\n\nfig : The figure or grid position object.\nres : One or more inv_out_1D structures containing the fit results.\n\nKeyword (optional) arguments:\n\ntitle : Title of the plot.\n\n\n\n\n\n","category":"method"},{"location":"functions/","page":"Functions","title":"Functions","text":"For expfits, we have:","category":"page"},{"location":"functions/","page":"Functions","title":"Functions","text":"plot(::NMRInversions.expfit_struct)\nplot!(::Union{Makie.Figure,Makie.GridPosition}, ::NMRInversions.expfit_struct )","category":"page"},{"location":"functions/#MakieCore.plot-Tuple{expfit_struct}","page":"Functions","title":"MakieCore.plot","text":"plot(res::expfit_struct...; kwargs...)\n\nPlot the results of an expfit call. \n\nThis function can take any number of expfit_struct structures as input. \n\ne.g. plot(data1) plots the data and the corresponding fit, but plot(data1, data2) or plot(data1, data2, data3) work as well, and they plot all of the results on the same plot.\n\nThis function creates a figure,  and calls plot!(fig, res) for each of the results in res, so please refer to the plot! documentation for more information.\n\n\n\n\n\n","category":"method"},{"location":"functions/#MakieCore.plot!-Tuple{Union{GridPosition, Figure}, expfit_struct}","page":"Functions","title":"MakieCore.plot!","text":"plot!(fig, res...; names, markersize, normeq)\n\nPlots the results of an expfit call on a figure or a grid position object.\n\nThe arguments are:\n\nfig : The figure or grid position object.\nres : One or more expfit_struct structures containing the fit results.\nnames : Vector with the names of the data (default is a vector of \"Data\" for each result).\nmarkersize : The size of the markers (default is 7).\nnormeq : Whether to plot the normalised form of the equation or not (default is true).\n\nNote that the res inputs are not a vector, but individual expfit_struct structures, like: plot!(fig, data1, data2, data3). If you want to use a vector of expfit_struct structures, make sure to  splat it by using ... in the function call (e.g. plot!(fig, [data1, data2, data3]...)).\n\n\n\n\n\n","category":"method"},{"location":"functions/","page":"Functions","title":"Functions","text":"And for 2D inversions, there's an interactive gui to characterize the inversion resutls:","category":"page"},{"location":"functions/","page":"Functions","title":"Functions","text":"plot(::NMRInversions.inv_out_2D)\nplot(::NMRInversions.inv_out_2D,::String)\nplot!(::Union{Makie.Figure,Makie.GridPosition}, ::NMRInversions.inv_out_2D )","category":"page"},{"location":"functions/#MakieCore.plot-Tuple{inv_out_2D}","page":"Functions","title":"MakieCore.plot","text":"Makie.plot(res::inv_out_2D)\n\nRun the GUI to plot the results and select peaks you want to label.\n\n\n\n\n\n","category":"method"},{"location":"functions/#MakieCore.plot-Tuple{inv_out_2D, String}","page":"Functions","title":"MakieCore.plot","text":"plot(results::NMRInversions.inv_out_2D, title::String; kwargs...)\n\nPlot the results contained in a inv_out_2D structure.\n\nThe arguments are:\n\nresults : The inv_out_2D structure containing the fit results.\ntitle : Title of the plot.\n\nkwargs are passed onto plot!.\n\n\n\n\n\n","category":"method"},{"location":"functions/#MakieCore.plot!-Tuple{Union{GridPosition, Figure}, inv_out_2D}","page":"Functions","title":"MakieCore.plot!","text":"plot!(fig, res::NMRInversions.inv_out_2D; title, clmap)\n\nPlot the results contained in a inv_out_2D structure on a figure or a grid position object.\n\nThe arguments are:\n\nfig : The figure or grid position object.\nres : The inv_out_2D structure containing the fit results.\n\nKeyword (optional) arguments:\n\ntitle : Title of the plot (default: \"\").\nclmap : Color map of the plot (default: :tempo).\n\n\n\n\n\n","category":"method"},{"location":"types_structs/#Types-and-Structures","page":"Types and Structures","title":"Types and Structures","text":"","category":"section"},{"location":"types_structs/","page":"Types and Structures","title":"Types and Structures","text":"In order to take advantage of julia's multiple dispatch,  we need to define some structures which can be used as input to our functions.","category":"page"},{"location":"types_structs/","page":"Types and Structures","title":"Types and Structures","text":"info: Info\nStructures are object which may contain multiple fields. The fields can be accessed using the dot syntax. For example, if we have a structure named foo, with two fields, a, and b,  we can access the value of a using foo.a and the value of b using foo.b.","category":"page"},{"location":"types_structs/#Pulse-sequences","page":"Types and Structures","title":"Pulse sequences","text":"","category":"section"},{"location":"types_structs/","page":"Types and Structures","title":"Types and Structures","text":"The following types are used so that the invert function  can create the appropriate kernel for the inversion.  Anywhere you see the seq keyword, one of the following must be used.","category":"page"},{"location":"types_structs/","page":"Types and Structures","title":"Types and Structures","text":"IR\nSR\nCPMG\nPFG\nIRCPMG","category":"page"},{"location":"types_structs/#NMRInversions.IR","page":"Types and Structures","title":"NMRInversions.IR","text":"Inversion recovery pulse sequence for 1D relaxation experiments. It can be used wherever the seq argument is required. \n\n\n\n\n\n","category":"type"},{"location":"types_structs/#NMRInversions.SR","page":"Types and Structures","title":"NMRInversions.SR","text":"Saturation recovery pulse sequence for 1D relaxation experiments. It can be used wherever the seq argument is required. \n\n\n\n\n\n","category":"type"},{"location":"types_structs/#NMRInversions.CPMG","page":"Types and Structures","title":"NMRInversions.CPMG","text":"CPMG pulse sequence for 1D relaxation experiments. It can be used wherever the seq argument is required. \n\n\n\n\n\n","category":"type"},{"location":"types_structs/#NMRInversions.PFG","page":"Types and Structures","title":"NMRInversions.PFG","text":"Pulsed field gradient pulse sequence for 1D diffusion experiments. It can be used wherever the seq argument is required. \n\n\n\n\n\n","category":"type"},{"location":"types_structs/#NMRInversions.IRCPMG","page":"Types and Structures","title":"NMRInversions.IRCPMG","text":"Inversion recovery - CPMG pulse sequence for 2D relaxation experiments (T1-T2). The direct dimension is the T2, and the indirect dimension is the T1 acquisition times. It can be used wherever the seq argument is required.\n\n\n\n\n\n","category":"type"},{"location":"types_structs/#Inversion-solvers","page":"Types and Structures","title":"Inversion solvers","text":"","category":"section"},{"location":"types_structs/","page":"Types and Structures","title":"Types and Structures","text":"These are used to let the invert fucntion know which solver to use. They can be used as input to the invert function as the 'solver' argument. (e.g., invert(data, solver=brd) or invert(data, solver=pdhgm(10,0.1) )).","category":"page"},{"location":"types_structs/","page":"Types and Structures","title":"Types and Structures","text":"brd\npdhgm\noptim_nnls","category":"page"},{"location":"types_structs/#NMRInversions.brd","page":"Types and Structures","title":"NMRInversions.brd","text":"brd\n\nSolver for tikhonov (L2) regularization, following this paper from Venkataramanan et al. Very fast, but only uses the identity as tiknohonov matrix. No options required, it just works. It can be used as a \"solver\" for the invert function.\n\n\n\n\n\n","category":"type"},{"location":"types_structs/#NMRInversions.pdhgm","page":"Types and Structures","title":"NMRInversions.pdhgm","text":"pdhgm(σ, τ)\n\nPrimal dual hybrid gradient method for L1 regularization,  following this paper from Reci et al. / Journal of Magnetic Resonance 281 (2017) 188–198 It can be used as a \"solver\" for the invert function.\n\nThe particular choice of σ and τ is heuristic.  A smaller σ will increase the stability while reducing the convergence speed of the algorithm. A good compromise between the two was found when σ = 0.1 and τ = 10.  The best values of σ and τ will depend slightly on the scaling of the signal.  Therefore, it is best to normalize the NMR signal to a maximum of 1, a technique which was followed in the cited study.\n\n\n\n\n\n","category":"type"},{"location":"types_structs/#NMRInversions.optim_nnls","page":"Types and Structures","title":"NMRInversions.optim_nnls","text":"optim_nnls(order)\n\nSimple non-negative least squares method for tikhonov (L2) regularization,  implemented using OptimizationOptimJl. All around effective, but can be slow for large problems, such as 2D inversions. It can be used as a \"solver\" for invert function. Order determines the tikhonov matrix. If 0 is chosen, the identity matrix is used.\n\n\n\n\n\n","category":"type"},{"location":"types_structs/#Finding-optimal-alpha","page":"Types and Structures","title":"Finding optimal alpha","text":"","category":"section"},{"location":"types_structs/","page":"Types and Structures","title":"Types and Structures","text":"These are methods for finding the optimal regularization parameter.  They can be used as input to the invert function as the 'alpha' argument (e.g., invert(data, alpha=gcv) ). If you'd like to use a particular value of alpha,  you can just use that number instead (invert(data, alpha=1).","category":"page"},{"location":"types_structs/","page":"Types and Structures","title":"Types and Structures","text":"gcv\nlcurve","category":"page"},{"location":"types_structs/#NMRInversions.gcv","page":"Types and Structures","title":"NMRInversions.gcv","text":"gcv\n\nGeneralized cross validation for finding the optimal regularization parameter α.\n\n\n\n\n\n","category":"type"},{"location":"types_structs/#NMRInversions.lcurve","page":"Types and Structures","title":"NMRInversions.lcurve","text":"lcurve\n\nL curve method for finding the optimal regularization parameter α.\n\nSTILL UNDER DEVELOPMENT!\n\n\n\n\n\n","category":"type"},{"location":"types_structs/#Inversion-and-expfit-inputs","page":"Types and Structures","title":"Inversion and expfit inputs","text":"","category":"section"},{"location":"types_structs/","page":"Types and Structures","title":"Types and Structures","text":"These can be used as inputs containing all of the necessary information to run the inversion and expfit functions  (e.g. if a is an input1D object, invert(a) does the job).","category":"page"},{"location":"types_structs/","page":"Types and Structures","title":"Types and Structures","text":"input1D\ninput2D","category":"page"},{"location":"types_structs/#NMRInversions.input1D","page":"Types and Structures","title":"NMRInversions.input1D","text":"input1D(seq, x, y)\n\nA structure containing the following fields:\n\nseq is the 1D pulse sequence (e.g. IR, CPMG, PFG)\nx, the x values of the measurement (e.g time for relaxation or b-factor for diffusion).\ny, the y values of the measurement. \n\nIt can be used as an input for the invert and expfit functions.\n\n\n\n\n\n","category":"type"},{"location":"types_structs/#NMRInversions.input2D","page":"Types and Structures","title":"NMRInversions.input2D","text":"input2D(seq, x, y)\n\nA structure containing the following fields:\n\nseq is the 2D pulse sequence (e.g. IRCPMG)\nx_direct is the direct dimension acquisition parameter (e.g. the times when you aquire CPMG echoes).\nx_indirect is the indirect dimension acquisition parameter (e.g. all the delay times τ in your IR sequence).\ndata is the 2D data matrix.\n\nIt can be used as an input for the invert function.\n\n\n\n\n\n","category":"type"},{"location":"types_structs/#Inversion-and-expfit-outputs","page":"Types and Structures","title":"Inversion and expfit outputs","text":"","category":"section"},{"location":"types_structs/","page":"Types and Structures","title":"Types and Structures","text":"When you run an inversion function (e.g. r = invert(a)), the output is an inv_out_1D or inv_out_2D structure.","category":"page"},{"location":"types_structs/","page":"Types and Structures","title":"Types and Structures","text":"inv_out_1D\ninv_out_2D\nexpfit_struct","category":"page"},{"location":"types_structs/#NMRInversions.inv_out_1D","page":"Types and Structures","title":"NMRInversions.inv_out_1D","text":"inv_out_1D(seq, x, y, xfit, yfit, X, f, r, SNR, α, wa)\n\nOutput of the invert function for 1D pulse sequences. A structure containing the following fields:\n\nseq is the 1D pulse sequence (e.g. IR, CPMG, PFG)\nx, the x values of the measurement (e.g time for relaxation or b-factor for diffusion).\ny, the y values of the measurement.\nxfit, the x values of the fitted data.\nyfit, the y values of the fitted data.\nX, the x values of the inversion results.\nf, the inversion results.\nr, the residuals.\nSNR, the signal-to-noise ratio.\nα, the regularization parameter.\nwa, the weighted average of the inversion results (e.g. the mean relaxation time or diffusion coefficient).\n\n\n\n\n\n","category":"type"},{"location":"types_structs/#NMRInversions.inv_out_2D","page":"Types and Structures","title":"NMRInversions.inv_out_2D","text":"inv_out_2D(seq, X_dir, X_indir, F, r, SNR, α, filter, selections)\n\nOutput of the invert function for 2D pulse sequences. A structure containing the following fields:\n\nseq is the 2D pulse sequence (e.g. IRCPMG)\nX_dir, the x values of the direct dimension.\nX_indir, the x values of the indirect dimension.\nF, the inversion results as a matrix.\nr, the residuals.\nSNR, the signal-to-noise ratio.\nalpha, the regularization parameter.\nfilter, apply a mask to filter out artefacts when plotting.\nselections, the selection masks (e.g. when you want to highlight some peaks in a T₁-T₂ map).\n\n\n\n\n\n","category":"type"},{"location":"types_structs/#NMRInversions.expfit_struct","page":"Types and Structures","title":"NMRInversions.expfit_struct","text":"Output of the expfit function. Structure containing information about multiexponential fits.\n\nThe fields are as follows:\n\nseq: The pulse sequence\nx : The x acquisition values (e.g. time for relaxation or b-factor for diffusion).\ny : The y acquisition values.\nu : The fitted parameters for the mexp function.\nu0 : The initial parameters for the mexp function.\nr : The residuals.\neq : The equation of the fitted function.\neqn : The equation of the initial function.\n\n\n\n\n\n","category":"type"},{"location":"types_structs/#Kernel-structure","page":"Types and Structures","title":"Kernel structure","text":"","category":"section"},{"location":"types_structs/","page":"Types and Structures","title":"Types and Structures","text":"svd_kernel_struct","category":"page"},{"location":"types_structs/#NMRInversions.svd_kernel_struct","page":"Types and Structures","title":"NMRInversions.svd_kernel_struct","text":"svd_kernel_struct(K,g,U,S,V)\n\nA structure containing the following fields:\n\nK, the kernel matrix.\nG, the data vector.\nU, the left singular values matrix.\nS, the singular values vector.\nV, the right singular values matrix.\n\nTo access the fields of a structure, we use the dot notation,  e.g. if the structure is called a and we want to access the kernel contained in there, we type a.K\n\n\n\n\n\n","category":"type"},{"location":"","page":"Overview","title":"Overview","text":"Welcome to the NMRInversions.jl package!","category":"page"},{"location":"","page":"Overview","title":"Overview","text":"Hopefully this documentation will provide all the information you need to get started.","category":"page"},{"location":"","page":"Overview","title":"Overview","text":"info: Installing julia\nYou can download the julia programming language by following this link. Afterwards, you can either use the julia console (aka REPL)  by finding the julia .exe file in your computer, or you can use VScode, which provides an environment similar to MATLAB with the use of the  julia extension (follow link for installation instructions).  Of course, if you already have a preferred development workflow  and you know what you're doing, by all means go for it.","category":"page"},{"location":"","page":"Overview","title":"Overview","text":"info: Installing the package\nThe package can be installed by running the following command on the julia console:using Pkg ; Pkg.add(\"NMRInversions\")This usually takes a while, but it needs to be done only once  (unless you swap environment, more on that  here).Afterwards, you can use the package by running using NMRInversionsin your julia console, every time you start a new session.Whenever a new version comes up, you can run:using Pkg ; Pkg.update(\"NMRInversions\")to update the package.","category":"page"},{"location":"","page":"Overview","title":"Overview","text":"info: GLMakie extension\nThe package provides an extension for interactive visualization, using the GLMakie package. To gain access to these capabilities, you also need to install GLMakie:using Pkg ; Pkg.add(\"GLMakie\")And to use it, you need to runusing GLMakiein order to access the plotting functions.","category":"page"},{"location":"","page":"Overview","title":"Overview","text":"For more details on how to use the package, you can start by refering to the tutorial section.","category":"page"},{"location":"tutorial/#Intro","page":"Tutorial","title":"Intro","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"This page will give you a basic idea about how to use the NMRInversions package. For more details, it's best to refer to the functions page.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"hint: Hint\nAll of the commands mentioned below should be typed in the julia console,  or saved in a text file with the .jl extension, to be used  from a terminal with julia file.jl, or through an IDE such as VSCode.","category":"page"},{"location":"tutorial/#Performing-an-inversion","page":"Tutorial","title":"Performing an inversion","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Suppose we're working with data coming from a Spinsolve instrument. Then we can do the following:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"using NMRInversions\n\ndata = import_spinsolve()","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"info: Info\nSince we called the import_spinsolve function without an argument,  it'll open a file dialog for us to select the files we want to import. (note that import_spinsolve requires two files, the aqcu.par file, plus the file containing the experiment data.)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Now we have the data imported, the inversion can be performed using a single line of code!","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"results = invert(data)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"info: Info\nThe results variable above is an inv_out_1D or inv_out_2D structure,  which contains all the relevant information produced by the inversion function. To access that information, we can look at the fields of the structure using the dot notation. The field names can be shown by using the REPL help mode (typing ? at the julia> prompt),  and typing the variable's name (in this case ?results).  Alternatively, running @doc results will also give you the same answers.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"The results can easily be visualised through the GLMakie extension of the package.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"using GLMakie\nplot(results)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"info: Info\nThe plot function of GLMakie is modified by this package  to work with results from the invert function as arguments. It's really easy to use, but if you want more control  on how your plots look, it's best to create them from scratch  using all the tools available in GLMakie.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Or, if the plot is all you need, one line of code is enough:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"using NMRInversions, GLMakie\nplot(invert(import_spinsolve()))","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Note that the workflow above can work for both 1D and 2D inversions!","category":"page"},{"location":"tutorial/#Using-the-expfit-function","page":"Tutorial","title":"Using the expfit function","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"In a similar way, we can perform various exponential  fits to the imported data using the expfit function.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"using NMRInversions, GLMakie\n\ndata = import_spinsolve()\n\na = expfit(1, data)  # mono-exponential fit\nb = expfit(2, data)  # bi-exponential fit\n\nplot(a,b)  # Visualize both on the same plot\n\na.eqn  # Print the equation of the mono-exponential fit\nb.eqn  # Print the equation of the bi-exponential fit","category":"page"}]
}
