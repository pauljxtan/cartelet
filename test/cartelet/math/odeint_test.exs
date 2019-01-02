defmodule OdeIntTest do
  use ExUnit.Case
  use Numerix.Tensor
  alias Cartelet.Math.OdeInt

  # Lorenz attractor test case
  @state_init Tensor.new([0.01, 0.01, 0.01])
  @ode &OdeInt.lorenz/3
  @dt 0.01
  @steps 10
  @lorenz_params [sigma: 10, beta: 8 / 3, rho: 28]

  # Manually computed results
  @expected %{t: 0.1, x: 0.021840, y: 0.046459, z: 0.007696}

  # Error tolerance for each method
  @euler_tol @dt / 2
  @heun_tol @dt / 50
  @midpoint_tol @dt / 50
  @ralston_tol @dt / 50
  @rk4_tol @dt / 10_000

  test "integrates Lorenz attractor with Euler's method" do
    {:ok, t, state} = OdeInt.integrate(@ode, @state_init, @dt, @steps)

    [x, y, z] = state.items
    assert_in_delta(t, @expected.t, @euler_tol)
    assert_in_delta(x, @expected.x, @euler_tol)
    assert_in_delta(y, @expected.y, @euler_tol)
    assert_in_delta(z, @expected.z, @euler_tol)

    # With optional arguments
    {:ok, t, state} =
      OdeInt.integrate(@ode, @state_init, @dt, @steps, :euler, @lorenz_params)

    [x, y, z] = state.items
    assert_in_delta(t, @expected.t, @euler_tol)
    assert_in_delta(x, @expected.x, @euler_tol)
    assert_in_delta(y, @expected.y, @euler_tol)
    assert_in_delta(z, @expected.z, @euler_tol)
  end

  test "integrates Lorenz attractor with Heun's method" do
    {:ok, t, state} = OdeInt.integrate(@ode, @state_init, @dt, @steps)

    [x, y, z] = state.items
    assert_in_delta(t, @expected.t, @heun_tol)
    assert_in_delta(x, @expected.x, @heun_tol)
    assert_in_delta(y, @expected.y, @heun_tol)
    assert_in_delta(z, @expected.z, @heun_tol)

    # With optional arguments
    {:ok, t, state} =
      OdeInt.integrate(@ode, @state_init, @dt, @steps, :heun, @lorenz_params)

    [x, y, z] = state.items
    assert_in_delta(t, @expected.t, @heun_tol)
    assert_in_delta(x, @expected.x, @heun_tol)
    assert_in_delta(y, @expected.y, @heun_tol)
    assert_in_delta(z, @expected.z, @heun_tol)
  end

  test "integrates Lorenz attractor with the midpoint method" do
    {:ok, t, state} = OdeInt.integrate(@ode, @state_init, @dt, @steps)

    [x, y, z] = state.items
    assert_in_delta(t, @expected.t, @midpoint_tol)
    assert_in_delta(x, @expected.x, @midpoint_tol)
    assert_in_delta(y, @expected.y, @midpoint_tol)
    assert_in_delta(z, @expected.z, @midpoint_tol)

    # With optional arguments
    {:ok, t, state} =
      OdeInt.integrate(
        @ode,
        @state_init,
        @dt,
        @steps,
        :midpoint,
        @lorenz_params
      )

    [x, y, z] = state.items
    assert_in_delta(t, @expected.t, @midpoint_tol)
    assert_in_delta(x, @expected.x, @midpoint_tol)
    assert_in_delta(y, @expected.y, @midpoint_tol)
    assert_in_delta(z, @expected.z, @midpoint_tol)
  end

  test "integrates Lorenz attractor with Raltson's method" do
    {:ok, t, state} = OdeInt.integrate(@ode, @state_init, @dt, @steps)

    [x, y, z] = state.items
    assert_in_delta(t, @expected.t, @ralston_tol)
    assert_in_delta(x, @expected.x, @ralston_tol)
    assert_in_delta(y, @expected.y, @ralston_tol)
    assert_in_delta(z, @expected.z, @ralston_tol)

    # With optional arguments
    {:ok, t, state} =
      OdeInt.integrate(
        @ode,
        @state_init,
        @dt,
        @steps,
        :ralston,
        @lorenz_params
      )

    [x, y, z] = state.items
    assert_in_delta(t, @expected.t, @ralston_tol)
    assert_in_delta(x, @expected.x, @ralston_tol)
    assert_in_delta(y, @expected.y, @ralston_tol)
    assert_in_delta(z, @expected.z, @ralston_tol)
  end

  test "integrates Lorenz attractor with RK4" do
    {:ok, t, state} = OdeInt.integrate(@ode, @state_init, @dt, @steps)

    [x, y, z] = state.items
    assert_in_delta(t, @expected.t, @rk4_tol)
    assert_in_delta(x, @expected.x, @rk4_tol)
    assert_in_delta(y, @expected.y, @rk4_tol)
    assert_in_delta(z, @expected.z, @rk4_tol)

    # With optional arguments
    {:ok, t, state} =
      OdeInt.integrate(@ode, @state_init, @dt, @steps, :rk4, @lorenz_params)

    [x, y, z] = state.items
    assert_in_delta(t, @expected.t, @rk4_tol)
    assert_in_delta(x, @expected.x, @rk4_tol)
    assert_in_delta(y, @expected.y, @rk4_tol)
    assert_in_delta(z, @expected.z, @rk4_tol)
  end
end
