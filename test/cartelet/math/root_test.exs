defmodule RootTest do
  use ExUnit.Case
  import Cartelet.Math.Root

  test "approximates the root of a univariate function" do
    f = fn x -> :math.pow(x - 1, 2) end
    root = 1.0
    # Perturbation constant
    delta = 1.0e-3
    target_err = 1.0e-3

    assert_in_delta(secant(f, 10, delta, target_err), root, target_err)

    # TODO: test edge cases
  end
end
