using Pkg

dependencies = [
  "ArgParse",
  "ColorSchemes",
  "DataFrames",
  "DelimitedFiles",
  "Gadfly",
  "StatsBase",
  "XAM",
]

Pkg.add(dependencies)

#One time precompilation of packages
create_sysimage([:XAM, :StatsBase, :DataFrames, :ColorSchemes, :Gadfly, :ArgParse, :DelimitedFiles], sysimage_path="m2img.so")
#Now we can start with julia --sysimage m2img.so
