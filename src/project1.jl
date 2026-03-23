#!/usr/bin/env julia
#
# Copyright 2022-2024 John T. Foster
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
module project1

using DifferentialEquations
using Plots
using LaTeXStrings

function lorenz_solver(x₀, ρ)
    σ = 10.0
    β = 8 / 3

    function lorenz!(du, u, p, t)
        x, y, z = u
        du[1] = σ * (y - x)
        du[2] = x * (ρ - z) - y
        du[3] = x * y - β * z
    end

    tspan = (0.0, 100.0)
    prob = ODEProblem(lorenz!, x₀, tspan)
    return solve(prob)
end

function lorenz_plot()
    x0_1 = [0.0, 1.0, 0.0]
    x0_2 = [1.0e-5, 1.0, 0.0]

    sol14_1 = lorenz_solver(x0_1, 14.0)
    sol14_2 = lorenz_solver(x0_2, 14.0)

    sol28_1 = lorenz_solver(x0_1, 28.0)
    sol28_2 = lorenz_solver(x0_2, 28.0)

    t14 = sol14_1.t
    t28 = sol28_1.t

    x14_1 = getindex.(sol14_1.u, 1)
    z14_1 = getindex.(sol14_1.u, 3)
    x14_2 = getindex.(sol14_2.u, 1)
    z14_2 = getindex.(sol14_2.u, 3)

    x28_1 = getindex.(sol28_1.u, 1)
    z28_1 = getindex.(sol28_1.u, 3)
    x28_2 = getindex.(sol28_2.u, 1)
    z28_2 = getindex.(sol28_2.u, 3)

    p1 = plot(
        t14, x14_1,
        title = L"\rho = 14.0",
        xlabel = "",
        ylabel = L"x_1, x_2",
        legend = false,
        grid = true
    )
    plot!(p1, t14, x14_2)

    p2 = plot(
        t28, x28_1,
        title = L"\rho = 28.0",
        xlabel = "",
        ylabel = "",
        legend = false,
        grid = true
    )
    plot!(p2, t28, x28_2)

    p3 = plot(
        t14, abs.(x14_1 .- x14_2),
        xlabel = L"t",
        ylabel = L"|x_1 - x_2|",
        legend = false,
        grid = true
    )

    p4 = plot(
        t28, abs.(x28_1 .- x28_2),
        xlabel = L"t",
        ylabel = "",
        legend = false,
        grid = true
    )

    p5 = plot(
        x14_1, z14_1,
        xlabel = L"x_1",
        ylabel = L"z_1",
        legend = false,
        grid = true
    )
    plot!(p5, x14_2, z14_2)

    p6 = plot(
        x28_1, z28_1,
        xlabel = L"x_1",
        ylabel = "",
        legend = false,
        grid = true
    )
    plot!(p6, x28_2, z28_2)

    return plot(
        p1, p2, p3, p4, p5, p6;
        layout = (3, 2),
        size = (900, 900),
        link = :x
    )
end

export lorenz_solver, lorenz_plot

end
