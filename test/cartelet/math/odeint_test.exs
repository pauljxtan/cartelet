defmodule OdeIntTest do
  use ExUnit.Case
  use Numerix.Tensor
  alias Cartelet.Math.OdeInt

  test "integrates Lorenz attractor with RK4" do
    state_init = Tensor.new([0.01, 0.01, 0.01])
    dt = 0.01
    steps = 10
    ode = &OdeInt.lorenz/3

    {:ok, t, state} = OdeInt.integrate(ode, state_init, dt, steps)

    # This is a pretty generous tolerance, but should be good enough for
    # testing program correctness
    tol = dt / 10

    [x, y, z] = state.items
    assert_in_delta(t, 0.1, tol)
    assert_in_delta(x, 0.021840, tol)
    assert_in_delta(y, 0.046459, tol)
    assert_in_delta(z, 0.007696, tol)

    # Specify optional arguments
    method = :rk4
    params = [sigma: 10, beta: 8 / 3, rho: 28]

    {:ok, t, state} =
      OdeInt.integrate(ode, state_init, dt, steps, method, params)

    [x, y, z] = state.items
    assert_in_delta(t, 0.1, tol)
    assert_in_delta(x, 0.021840, tol)
    assert_in_delta(y, 0.046459, tol)
    assert_in_delta(z, 0.007696, tol)
  end
end
