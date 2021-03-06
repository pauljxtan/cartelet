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
    k1 = dt * ode.(t, state, params)
    k2 = dt * ode.(t + dt, state + k1, params)

    t_new = t + dt
    state_new = state + k1 / 2 + k2 / 2
    {t_new, state_new}
  end
end
