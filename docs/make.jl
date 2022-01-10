push!(LOAD_PATH,joinpath(@__DIR__, ".."))
using Documenter, VacBEM

makedocs(
    modules = [VacBEM],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "blaise-faugeras",
    sitename = "VacBEM.jl",
    pages = Any["index.md"]
    # strict = true,
    # clean = true,
    # checkdocs = :exports,
)

deploydocs(
    repo = "github.com/blaise-faugeras/VacBEM.jl.git",
)
