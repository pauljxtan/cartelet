defmodule IntegralTest do
  use ExUnit.Case
  import Cartelet.Math.Integral

  test "evaluates the integral of a univariate function" do
    f = fn x -> 2 * :math.pow(x, 3) end
    fint = fn x -> 0.5 * :math.pow(x, 4) end
    a = 1
    b = 2
    expected = fint.(b) - fint.(a)

    # Trapezoidal rules
    assert_in_delta(trapezoid(f, a, b), expected, 1.0)
    assert_in_delta(trapezoid_composite(f, a, b, 1000), expected, 0.1)
  end
end
