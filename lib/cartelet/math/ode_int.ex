defmodule Cartelet.Math.OdeInt do
  @moduledoc """
  Provides functions for numerical integration of systems of ordinary
  differential equations.
  """
  use Numerix.Tensor
  @type state :: %Tensor{}
  @type ode_func :: (float, state, map -> state)

  # Start with just RK4, then generalize to other methods via defined behaviour
  @spec integrate(state, float, integer, ode_func, map) :: state
  def integrate(state_init, dt, steps, ode_func, ode_kwargs \\ %{}) do
    1..steps
    |> Enum.reduce({0, state_init}, fn _, {t, state} ->
      rk4(t, state, dt, ode_func, ode_kwargs)
    end)
  end

  # https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods
  defp rk4(t, state, dt, ode, kwargs) do
    k1 = dt * ode.(t, state, kwargs)
    k2 = dt * ode.(t + dt / 2, state + k1 / 2, kwargs)
    k3 = dt * ode.(t + dt / 2, state + k2 / 2, kwargs)
    k4 = dt * ode.(t + dt, state + k3, kwargs)

    t_new = t + dt
    state_new = state + (k1 + 2 * k2 + 2 * k3 + k4) / 6
    {t_new, state_new}
  end

  @lorenz_kwargs_default %{sigma: 10, beta: 8 / 3, rho: 28}
  @doc """
  Defines the Lorenz system (a.k.a. Lorenz attractor).
  https://en.wikipedia.org/wiki/Lorenz_system
  """
  def lorenz(_t, state, kwargs) do
    kwargs = if %{}, do: @lorenz_kwargs_default, else: kwargs
    %{sigma: sigma, beta: beta, rho: rho} = kwargs

    [x, y, z] = state.items
    x_new = sigma * (y - x)
    y_new = x * (rho - z) - y
    z_new = x * y - beta * z
    Tensor.new([x_new, y_new, z_new])
  end
end
