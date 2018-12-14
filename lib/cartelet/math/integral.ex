defmodule Cartelet.Math.Integral do
  @moduledoc """
  Provides functions for evaluating integrals of univariate functions.
  """
  require Integer

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
    {h, xs} = increments(a, b, n)
    first = xs |> Enum.at(0) |> f.()
    mid = xs |> Enum.slice(1, n - 1) |> Enum.map(f) |> Enum.sum()
    last = xs |> Enum.at(n - 1) |> f.()
    h / 2 * (first + 2 * mid + last)
  end

  @doc """
  Integrates a univariate `f` between `a` and `b` using Simpson's 1/3 rule.
  """
  def simpson_13(f, a, b) do
    {h, xs} = increments(a, b, 2)
    [x0, x1, x2] = xs
    h / 3 * (f.(x0) + 4 * f.(x1) + f.(x2))
  end

  @doc """
  Integrates a univariate `f` between `a` and `b` using Simpson's composite
  1/3 rule with `n` segments.
  """
  def simpson_13_composite(_f, _a, _b, n) when not Integer.is_even(n),
    do: "n must be even"

  def simpson_13_composite(f, a, b, n) do
    {h, xs} = increments(a, b, n)
    first = xs |> Enum.at(0) |> f.()
    odd_indices = 1..(n - 1) |> Enum.take_every(2)
    even_indices = 2..(n - 1) |> Enum.take_every(2)

    odd =
      odd_indices
      |> Enum.map(fn i -> Enum.at(xs, i) end)
      |> Enum.map(f)
      |> Enum.sum()

    even =
      even_indices
      |> Enum.map(fn i -> Enum.at(xs, i) end)
      |> Enum.map(f)
      |> Enum.sum()

    last = xs |> Enum.at(n - 1) |> f.()
    h / 3 * (first + 4 * odd + 2 * even + last)
  end

  @doc """
  Integrates a univariate `f` between `a` and `b` using Simpson's 3/8 rule.
  """
  def simpson_38(f, a, b) do
    {h, xs} = increments(a, b, 3)
    [x0, x1, x2, x3] = xs
    h * (3 / 8) * (f.(x0) + 3 * f.(x1) + 3 * f.(x2) + f.(x3))
  end

  defp increments(a, b, n) do
    h = (b - a) / n
    xs = 0..n |> Enum.map(fn i -> a + i * h end)
    {h, xs}
  end
end
