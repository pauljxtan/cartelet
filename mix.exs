defmodule Cartelet.MixProject do
  use Mix.Project

  def project do
    [
      app: :cartelet,
      version: "0.1.0",
      elixir: "~> 1.7",
      start_permanent: Mix.env() == :prod,
      deps: deps()
    ]
  end

  # Run "mix help compile.app" to learn about applications.
  def application do
    [
      extra_applications: [:logger]
    ]
  end

  # Run "mix help deps" to learn about dependencies.
  defp deps do
    [
      # {:dep_from_hexpm, "~> 0.3.0"},
      # {:dep_from_git, git: "https://github.com/elixir-lang/my_dep.git", tag: "0.1.0"},
      {:numerix, "~> 0.5"},
      {:ex_doc, "~> 0.19", only: :dev, runtime: false},
      # {:dialyxir, "~> 0.5", only: [:dev], runtime: false}
      {:dialyxir, "~> 1.0.0-rc.4", only: [:dev], runtime: false}
    ]
  end
end
