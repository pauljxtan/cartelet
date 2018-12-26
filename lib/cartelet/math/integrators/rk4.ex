defmodule Cartelet.Math.Integrators.Rk4 do
  @moduledoc """
  Implements the ODE integrator behaviour for 4th-order Runge-Kutta (RK4).
  """
  use Numerix.Tensor

  @behaviour Cartelet.Math.Integrators.OdeIntegrator

  @doc """
  Integrates an ODE system over a single step, at a given point in time, via
  the 4th-order Runge-Kutta method (RK4).

  Reference: [Wikipedia](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods)
  """
  @impl Cartelet.Math.Integrators.OdeIntegrator
  @spec step(OdeInt.ode_func(), float, OdeInt.state(), float, keyword) ::
          {float, OdeInt.state()}
  def step(ode, t, state, dt, params) do
    k1 = dt * ode.(t, state, params)
    k2 = dt * ode.(t + dt / 2, state + k1 / 2, params)
    k3 = dt * ode.(t + dt / 2, state + k2 / 2, params)
    k4 = dt * ode.(t + dt, state + k3, params)

    t_new = t + dt
    state_new = state + (k1 + 2 * k2 + 2 * k3 + k4) / 6
    {t_new, state_new}
  end
end
