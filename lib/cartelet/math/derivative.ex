defmodule Cartelet.Math.Derivative do
  @moduledoc """
  Provides functions evaluating derivatives of univariate functions.

  https://en.wikipedia.org/wiki/Finite_difference
  """
  # TODO: Include error estimates for each method?
  # TODO: Generalized nth-order derivatives

  @doc """
  Evaluates the first derivative of a univariate function `f` at the value `x`,
  with step size `h`, using a forward difference scheme.
  """
  @spec forward_difference((float -> float), float, float) :: float
  def forward_difference(f, x, h), do: (f.(x + h) - f.(x)) / h

  @doc """
  Evaluates the first derivative of a univariate function `f` at the value `x`,
  with step size `h`, using a backward difference scheme.
  """
  @spec backward_difference((float -> float), float, float) :: float
  def backward_difference(f, x, h), do: (f.(x) - f.(x - h)) / h

  @doc """
  Evaluates the first derivative of a univariate function `f` at the value `x`,
  with step size `h`, using a central difference scheme.
  """
  @spec central_difference((float -> float), float, float) :: float
  def central_difference(f, x, h), do: (f.(x + h / 2) - f.(x - h / 2)) / h

  @doc """
  Evaluates the second derivative of a univariate function `f` at the value `x`,
  with step size `h`, using a forward difference scheme.
  """
  @spec forward_difference_second((float -> float), float, float) :: float
  def forward_difference_second(f, x, h),
    do: (f.(x + 2 * h) - 2 * f.(x + h) + f.(x)) / :math.pow(h, 2)

  @doc """
  Evaluates the second derivative of a univariate function `f` at the value `x`,
  with step size `h`, using a backward difference scheme.
  """
  @spec backward_difference_second((float -> float), float, float) :: float
  def backward_difference_second(f, x, h),
    do: (f.(x) - 2 * f.(x - h) + f.(x - 2 * h)) / :math.pow(h, 2)

  @doc """
  Evaluates the second derivative of a univariate function `f` at the value `x`,
  with step size `h`, using a central difference scheme.
  """
  @spec central_difference_second((float -> float), float, float) :: float
  def central_difference_second(f, x, h),
    do: (f.(x + h) - 2 * f.(x) + f.(x - h)) / :math.pow(h, 2)
end
