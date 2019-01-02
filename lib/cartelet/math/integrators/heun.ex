defmodule Cartelet.Math.Integrators.Heun do
  @moduledoc """
  Implements the ODE integrator behaviour for Heun's method.
  """
  use Numerix.Tensor

  @behaviour Cartelet.Math.Integrators.OdeIntegrator

  @doc """
  Integrates an ODE system over a single step, at a given point in time, via
  Heun's method.

  Reference: [Wikipedia](https://en.wikipedia.org/wiki/Heun%27s_method)
  """
  @impl Cartelet.Math.Integrators.OdeIntegrator
  @spec step(OdeInt.ode_func(), float, OdeInt.state(), float, keyword) ::
          {float, OdeInt.state()}
  def step(ode, t, state, dt, params) do
    t_new = t + dt
    predictor = state + dt * ode.(t, state, params)
    avg_slope = (ode.(t, state, params) + ode.(t_new, predictor, params)) / 2
    state_new = state + dt * avg_slope
    {t_new, state_new}
  end
end
