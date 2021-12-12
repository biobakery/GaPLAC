# ./gaplac -v sample "y :~| SqExp(:x; l=1.5)" --at "x = rand(Uniform(-5,5), 50)" --output data_sqexp.tsv --plot notes/assets/sqexpplot.png

using GaPLAC
GaPLAC.runtests("Interface")

using GaPLAC.AbstractGPs
using GaPLAC.CairoMakie

cli_args = Dict(
    "spec"   => "y :~| SqExp(:x; l=1.5)",
    "at"     => "x = rand(Uniform(-5,5), 50)",
    "output" => "data_sqexp.tsv",
    "plot"   => "notes/assets/sqexpplot.png"
)

gpspec = GaPLAC.gp_spec(cli_args["spec"])
gp, vars = GaPLAC.make_gp(gpspec)

atdict = GaPLAC.getatrange(cli_args["at"], vars)


df = GaPLAC._make_test_df((atdict[v] for v in vars)...; vars)
X = RowVecs(Matrix(df))
df[!, GaPLAC.response(gpspec)] = rand(gp(X, 0.1))

GaPLAC._df_output(df, cli_args)

fig, ax, l = GaPLAC.sample_plot(gpspec, df)
fig
save(cli_args["plot"], fig)

# ./gaplac -v sample "y :~| OU(:x; l=1.5)" --at "x = rand(Uniform(-5,5), 50)" --output data_ou.tsv --plot notes/assets/ouplot.png

cli_args = Dict(
    "spec"   => "y :~| OU(:x; l=1.5)",
    "at"     => "x = rand(Uniform(-5,5), 50)",
    "output" => "data_ou.tsv",
    "plot"   => "notes/assets/ouplot.png"
)

gpspec = GaPLAC.gp_spec(cli_args["spec"])
gp, vars = GaPLAC.make_gp(gpspec)

atdict = GaPLAC.getatrange(cli_args["at"], vars)


df = GaPLAC._make_test_df((atdict[v] for v in vars)...; vars)
X = RowVecs(Matrix(df))
df[!, GaPLAC.response(gpspec)] = rand(gp(X, 0.1))

GaPLAC._df_output(df, cli_args)

fig, ax, l = GaPLAC.sample_plot(gpspec, df)
fig
save(cli_args["plot"], fig)

# ./gaplac -v sample "y :~| SqExp(:x)" --at "x = rand(Uniform(-5,5), 50)" --output data_ou.tsv --plot notes/assets/ouplot.png