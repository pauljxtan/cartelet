defmodule Cartelet.Math.Integral do
  @moduledoc """
  Provides functions for evaluating integrals of univariate functions.
  """
  @doc """
  Integrates a univariate function `f` between `a` and `b` using the
  trapezoidal rule.
  """
  @spec trapezoid((float -> float), float, float) :: float
  def trapezoid(f, a, b), do: (b - a) / 2 * (f.(b) - f.(a))

  @doc """
  Integrates a univariate function `f` between `a` and `b` using the
  composite trapezoidal rule with `n` segments.
  """
  @spec trapezoid_composite((float -> float), float, float, integer) :: float
  def trapezoid_composite(f, a, b, n \\ 1000) do
    h = (b - a) / n
    xs = 1..n |> Enum.map(fn i -> a + i * h end)
    first = xs |> Enum.at(0) |> f.()
    last = xs |> Enum.at(n - 1) |> f.()
    mid = xs |> Enum.slice(1, n - 1) |> Enum.map(f) |> Enum.sum()
    h / 2 * (first + 2 * mid + last)
  end
end
