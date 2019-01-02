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

  @typedoc """
  Represents an integration method.
  """
  @type method :: :euler | :heun | :midpoint | :ralston | :rk4

  @integrator_lookup [
    {:euler, Cartelet.Math.Integrators.Euler},
    {:heun, Cartelet.Math.Integrators.Heun},
    {:midpoint, Cartelet.Math.Integrators.Midpoint},
    {:ralston, Cartelet.Math.Integrators.Ralston},
    {:rk4, Cartelet.Math.Integrators.Rk4}
  ]

  @doc """
  Integrates the given ODE system with the given initial conditions and step
  size, over a number of steps.

  If the ODE parameters and integration method are not specified, defaults to
  pre-defined parameters and 4th-order Runge-Kutta.

  Returns the new state.
  """
  @spec integrate(ode_func, state, float, integer, method, keyword) ::
          {:ok, float, state}
  def integrate(ode, state_init, dt, steps, method \\ :rk4, params \\ []) do
    {t, state} = integrate_do(ode, 0.0, state_init, dt, steps, method, params)
    {:ok, t, state}
  end

  @spec integrate_do(ode_func, float, state, float, integer, method, keyword) ::
          {float, state}
  defp integrate_do(_, t, state, _, 0, _, _), do: {t, state}

  defp integrate_do(ode, t, state, dt, steps, method, params) do
    {t, state} = @integrator_lookup[method].step(ode, t, state, dt, params)
    integrate_do(ode, t, state, dt, steps - 1, method, params)
  end

  @doc """
  Evaluates the Lorenz system (a.k.a. Lorenz attractor).

  Returns the new state.

  Reference: [Wikipedia](https://en.wikipedia.org/wiki/Lorenz_system)
  """
  @lorenz_params_default [sigma: 10, beta: 8 / 3, rho: 28]
  @spec lorenz(float, state, keyword) :: state
  def lorenz(_t, state, params) do
    params = if params == [], do: @lorenz_params_default, else: params

    [x, y, z] = state.items
    x_new = params[:sigma] * (y - x)
    y_new = x * (params[:rho] - z) - y
    z_new = x * y - params[:beta] * z
    Tensor.new([x_new, y_new, z_new])
  end
end
