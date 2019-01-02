defmodule Cartelet.Math.Integrators.Euler do
  @moduledoc """
  Implements the ODE integrator behaviour for Euler's method.
  """
  use Numerix.Tensor

  @behaviour Cartelet.Math.Integrators.OdeIntegrator

  @doc """
  Integrates an ODE system over a single step, at a given point in time, via
  Euler's method.

  Reference: [Wikipedia](https://en.wikipedia.org/wiki/Euler_method)
  """
  @impl Cartelet.Math.Integrators.OdeIntegrator
  @spec step(OdeInt.ode_func(), float, OdeInt.state(), float, keyword) ::
          {float, OdeInt.state()}
  def step(ode, t, state, dt, params) do
    t_new = t + dt
    state_new = state + dt * ode.(t, state, params)
    {t_new, state_new}
  end
end
