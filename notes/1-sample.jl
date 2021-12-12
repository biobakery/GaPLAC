# ./gaplac -v sample "y :~| SqExp(:t; l=1.5)" --at "t = rand(Uniform(1,10), 50)" --output data_sqexp.tsv --plot notes/assets/sqexpplot.png

using GaPLAC
GaPLAC.runtests("Interface")

using GaPLAC.AbstractGPs
using GaPLAC.CairoMakie

cli_args = Dict(
    "spec"   => "y :~| SqExp(:t; l=1.5)",
    "at"     => "t = rand(Uniform(1,10), 50)",
    "output" => "notes/assets/data_sqexp.tsv",
    "plot"   => "notes/assets/sqexpplot.png"
)

gpspec = GaPLAC.gp_spec(cli_args["spec"])
gp, vars = GaPLAC.make_gp(gpspec)

atdict = GaPLAC.getatrange(gpspec, cli_args["at"], vars)


df = GaPLAC._make_test_df((atdict[v] for v in vars)...; vars)
X = RowVecs(Matrix(df))
df[!, GaPLAC.response(gpspec)] = rand(gp(X, 0.1))

GaPLAC._df_output(df, cli_args)

fig, ax, l = GaPLAC.sample_plot(gpspec, df)
fig
save(cli_args["plot"], fig)

# ./gaplac -v sample "y :~| OU(:t; l=1.5)" --at "t = rand(Uniform(1,10), 50)" --output data_ou.tsv --plot notes/assets/ouplot.png

cli_args = Dict(
    "spec"   => "y :~| OU(:t; l=1.5)",
    "at"     => "t = rand(Uniform(1,10), 50)",
    "output" => "notes/assets/data_ou.tsv",
    "plot"   => "notes/assets/ouplot.png"
)

gpspec = GaPLAC.gp_spec(cli_args["spec"])
gp, vars = GaPLAC.make_gp(gpspec)

atdict = GaPLAC.getatrange(gpspec, cli_args["at"], vars)


df = GaPLAC._make_test_df((atdict[v] for v in vars)...; vars)
X = RowVecs(Matrix(df))
df[!, GaPLAC.response(gpspec)] = rand(gp(X, 0.1))

GaPLAC._df_output(df, cli_args)

fig, ax, l = GaPLAC.sample_plot(gpspec, df)
fig
save(cli_args["plot"], fig)

# ./gaplac -v sample "y :~| SqExp(:t) + Linear(:x)" --at "t = rand(Uniform(1,10), 50)" --output data_multi.tsv --plot notes/assets/multiplot.png

cli_args = Dict(
    "spec"   => "y :~| SqExp(:t) + Linear(:x)",
    "at"     => "t = rand(Uniform(1,10), 50)",
    "output" => "notes/assets/data_multi.tsv",
    "plot"   => "notes/assets/multiplot.png"
)

gpspec = GaPLAC.gp_spec(cli_args["spec"])
gp, vars = GaPLAC.make_gp(gpspec)

atdict = GaPLAC.getatrange(gpspec, cli_args["at"], vars)

df = GaPLAC._make_test_df((atdict[v] for v in vars)...; vars)
X = RowVecs(Matrix(df))
df[!, GaPLAC.response(gpspec)] = rand(gp(X, 0.1))

GaPLAC._df_output(df, cli_args)

fig, ax, l = GaPLAC.sample_plot(gpspec, df)
fig
save(cli_args["plot"], fig)