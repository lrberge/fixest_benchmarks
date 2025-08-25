using DataFrames
using StatsModels
using FixedEffectModels
using GLFixedEffectModels
using Vcov

function jl_feols_timer(data::DataFrame, fml::String; vcov::String="", nthreads::Int=2)
  start_time = time()

  # Convert string to formula if needed
  formula = eval(Meta.parse("@formula(" * fml * ")"))
  if vcov === ""
    _ = reg(data, formula, nthreads=nthreads)
  else
    _ = reg(data, formula, Vcov.cluster(Base.Symbol(vcov)), nthreads=nthreads)
  end

  elapsed_time = time() - start_time
  return elapsed_time
end

function jl_poisson_timer(data::DataFrame, fml::String; vcov=nothing)
  start_time = time()

  # Convert string to formula if needed
  formula = eval(Meta.parse("@formula(" * fml * ")"))
  if vcov === nothing
    _ = nlreg(data, formula, Poisson(), LogLink())
  else
    _ = nlreg(data, formula, Poisson(), LogLink(), Vcov.cluster(:fe1))
  end

  elapsed_time = time() - start_time
  return elapsed_time
end

function jl_logit_timer(data::DataFrame, fml::String; vcov=nothing)
  start_time = time()

  # Convert string to formula if needed
  formula = eval(Meta.parse("@formula(" * fml * ")"))
  if vcov === nothing
    _ = nlreg(data, formula, Binomial(), LogitLink())
  else
    _ = nlreg(data, formula, Binomial(), LogitLink(), Vcov.cluster(:fe1))
  end

  elapsed_time = time() - start_time
  return elapsed_time
end




