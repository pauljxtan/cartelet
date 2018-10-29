defmodule CarteletTest do
  use ExUnit.Case
  doctest Cartelet

  test "greets the world" do
    assert Cartelet.hello() == :world
  end
end
