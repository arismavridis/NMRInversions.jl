module test_ext

using NMRInversions, CSV

function NMRInversions.foo(a::NMRInversions.ip_solver) 
    display("solver is ip_solver")
end

export foo

end # module