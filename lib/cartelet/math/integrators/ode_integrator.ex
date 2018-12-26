defmodule Cartelet.Math.Integrators.OdeIntegrator do
  @moduledoc """
  Provides a behaviour for ODE integrators.
  """
  @callback step(OdeInt.ode_func(), float, OdeInt.state(), float, keyword) ::
              {float, OdeInt.state()}
end
