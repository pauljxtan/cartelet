defmodule OdeIntTest do
  use ExUnit.Case

  use Numerix.Tensor
  alias Cartelet.Math.OdeInt

  test "integrates lorenz attractor" do
    state = Tensor.new([0.01, 0.01, 0.01])
    dt = 0.01
    steps = 10
    ode = &OdeInt.lorenz/3

    {t, state} = OdeInt.integrate(state, dt, steps, ode)

    # This is a pretty generous tolerance, but should be good enough for
    # testing program correctness
    tol = dt / 10

    [x, y, z] = state.items
    assert_in_delta(t, 0.1, tol)
    assert_in_delta(x, 0.021840, tol)
    assert_in_delta(y, 0.046459, tol)
    assert_in_delta(z, 0.007696, tol)
  end
end
