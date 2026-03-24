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

    u0 = [x₀, 1.0, 0.0]
    tspan = (0.0, 100.0)

    function lorenz(u, p, t)
        x, y, z = u
        return [
            σ * (y - x),
            x * (ρ - z) - y,
            x * y - β * z
        ]
    end

    prob = ODEProblem(lorenz, u0, tspan)
    return solve(prob)
end

function lorenz_plot()
    sol14_1 = lorenz_solver(0.0, 14)
    sol14_2 = lorenz_solver(1.0e-5, 14)

    sol28_1 = lorenz_solver(0.0, 28)
    sol28_2 = lorenz_solver(1.0e-5, 28)

    x14_1 = getindex.(sol14_1.u, 1)
    z14_1 = getindex.(sol14_1.u, 3)
    x14_2 = getindex.(sol14_2.u, 1)

    x28_1 = getindex.(sol28_1.u, 1)
    z28_1 = getindex.(sol28_1.u, 3)
    x28_2 = getindex.(sol28_2.u, 1)

    p1 = plot(
        sol14_1.t, x14_1,
        title = L"\rho = 14.0",
        ylabel = L"x_1,\ x_2",
        legend = false
    )
    plot!(p1, sol14_2.t, x14_2)

    p2 = plot(
        sol28_1.t, x28_1,
        title = L"\rho = 28.0",
        ylabel = "",
        legend = false
    )
    plot!(p2, sol28_2.t, x28_2)

    p3 = plot(
        sol14_1.t, abs.(x14_1 .- x14_2),
        xlabel = L"t",
        ylabel = L"|x_1 - x_2|",
        legend = false
    )

    p4 = plot(
        sol28_1.t, abs.(x28_1 .- x28_2),
        xlabel = L"t",
        ylabel = "",
        legend = false
    )

    p5 = plot(
        x14_1, z14_1,
        xlabel = L"x_1",
        ylabel = L"z_1",
        legend = false
    )

    p6 = plot(
        x28_1, z28_1,
        xlabel = L"x_1",
        ylabel = "",
        legend = false
    )

    return plot(
        p1, p2, p3, p4, p5, p6;
        layout = (3, 2),
        link = :x
    )
end

export lorenz_solver, lorenz_plot

end
