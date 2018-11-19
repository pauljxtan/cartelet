defmodule Cartelet.Math.Root do
  @moduledoc """
  Provides functions for finding roots of functions.
  """
  def secant(f, guess, delta, target_err) do
    x =
      guess -
        delta * guess * f.(guess) / (f.(guess + delta * guess) - f.(guess))

    err = abs((x - guess) / x)
    if err <= target_err, do: x, else: secant(f, x, delta, target_err)
  end
end
