defmodule DerivativeTest do
  use ExUnit.Case
  import Cartelet.Math.Derivative

  test "evaluates the derivatives of a univariate function" do
    f = fn x -> 2 * :math.pow(x, 3) end
    fp = fn x -> 6 * :math.pow(x, 2) end
    fpp = fn x -> 12 * x end
    x = 4

    # First-order derivative
    assert_in_delta(forward_difference(f, x, 1.0e-6), fp.(x), 1.0e-3)
    assert_in_delta(backward_difference(f, x, 1.0e-6), fp.(x), 1.0e-3)
    assert_in_delta(central_difference(f, x, 1.0e-6), fp.(x), 1.0e-3)

    # Second-order derivative
    assert_in_delta(forward_difference_second(f, x, 1.0e-6), fpp.(x), 1.0e-1)
    assert_in_delta(backward_difference_second(f, x, 1.0e-6), fpp.(x), 1.0e-1)
    assert_in_delta(central_difference_second(f, x, 1.0e-6), fpp.(x), 1.0e-2)
  end
end
