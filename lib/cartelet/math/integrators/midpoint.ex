defmodule Cartelet.Math.Integrators.Midpoint do
  @moduledoc """
  Implements the ODE integrator behaviour for the midpoint method.
  """
  use Numerix.Tensor

  @behaviour Cartelet.Math.Integrators.OdeIntegrator

  @doc """
  Integrates an ODE system over a single step, at a given point in time, via
  the midpoint method.

  Reference: [Wikipedia](https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods#Explicit_midpoint_method)
  """
  @impl Cartelet.Math.Integrators.OdeIntegrator
  @spec step(OdeInt.ode_func(), float, OdeInt.state(), float, keyword) ::
          {float, OdeInt.state()}
  def step(ode, t, state, dt, params) do
    k1 = dt * ode.(t, state, params)
    k2 = dt * ode.(t + dt / 2, state + k1 / 2, params)

    t_new = t + dt
    state_new = state + k2
    {t_new, state_new}
  end
end
