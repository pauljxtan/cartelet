defmodule Cartelet.Math.OdeInt do
  @moduledoc """
  Provides functions for numerical integration of systems of ordinary
  differential equations.

  (N.B.: "System" in this context can refer to either systems of ODEs or the
  system itself that we are integrating. We'll stick with this slightly
  confusing terminology until I come up with something better...)
  """
  use Numerix.Tensor

  @typedoc """
  Represents the state of a system of any dimension.
  """
  @type state :: %Tensor{}
  @typedoc """
  Evaluates a system of ODEs on a state at a given point in time.
  """
  @type ode_func :: (float, state, keyword -> state)

  @doc """
  Integrates the given ODE system with the given initial conditions and step
  size, over a number of steps.

  If the ODE parameters and integration method are not specified, defaults to
  pre-defined parameters and 4th-order Runge-Kutta.

  Returns the new state.
  """
  @spec integrate(ode_func, state, float, integer, keyword, atom) :: state
  def integrate(ode, state_init, dt, steps, params \\ [], method \\ :rk4) do
    {t, state} =
      1..steps
      |> Enum.reduce({0, state_init}, fn _, {t, state} ->
        case method do
          :rk4 -> rk4(ode, t, state, dt, params)
          _ -> :method_not_implemented
        end
      end)

    {:ok, t, state}
  end

  @doc """
  Integrates an ODE system over a single step, at a given point in time, via
  the 4th-order Runge-Kutta method (RK4).

  Reference: [Wikipedia](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods)
  """
  @spec rk4(ode_func, float, state, float, keyword) :: {float, state}
  def rk4(ode, t, state, dt, params) do
    k1 = dt * ode.(t, state, params)
    k2 = dt * ode.(t + dt / 2, state + k1 / 2, params)
    k3 = dt * ode.(t + dt / 2, state + k2 / 2, params)
    k4 = dt * ode.(t + dt, state + k3, params)

    t_new = t + dt
    state_new = state + (k1 + 2 * k2 + 2 * k3 + k4) / 6
    {t_new, state_new}
  end

  # TODO: Add more methods

  @doc """
  Evaluates the Lorenz system (a.k.a. Lorenz attractor).

  Returns the new state.

  Reference: [Wikipedia](https://en.wikipedia.org/wiki/Lorenz_system)
  """
  @lorenz_params_default [sigma: 10, beta: 8 / 3, rho: 28]
  @spec lorenz(float, state, keyword) :: state
  def lorenz(_t, state, params) do
    params = if %{}, do: @lorenz_params_default, else: params

    [x, y, z] = state.items
    x_new = params[:sigma] * (y - x)
    y_new = x * (params[:rho] - z) - y
    z_new = x * y - params[:beta] * z
    Tensor.new([x_new, y_new, z_new])
  end
end
