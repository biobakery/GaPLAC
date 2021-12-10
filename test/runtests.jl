using ReTest
include("GaPLACTests.jl")
GaPLACTests.runtests()

using GaPLAC
GaPLAC.runtests()