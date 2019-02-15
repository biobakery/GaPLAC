
using Distributions

abstract type AbstractGP
end

struct LaplaceGP <: AbstractGP
    # The covariance function
	# Signature: (x_i::Vector, x_j::Vector, self::Boolean, θc::Vector) -> Real
    cf::Function # Covariance function used for prediction
	cftr::Function # Covariance function used for training
	# Covariance function hyperparameter transformation
    θc_link::Function
    # Prior on covariance function hyperparameters
    θc_prior::Vector{UnivariateDistribution}
	# Covariance function hyperparameter names
	θc_names::Vector{String}
    # Jitter
    jitter::Float64

    # Data Likelihood
	# Signature: (f_i::Float64, z_i::Vector, θl::Vector) -> UnivariateDistribution
	datalik::Function
    # Data likelihood hyperparameter transformation
    θl_link::Function
    # Prior on data likelihood hyperparameters
    θl_prior::Vector{UnivariateDistribution}
	# Data likelihood hyperparameter names
	θl_names::Vector{String}
end
