#!/usr/bin/env julia
#
# Copyright 2022-2024 John T. Foster
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
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

tickfmt(x) = string(round(x, digits = 1))
blankfmt(x) = ""

function lorenz_solver(x₀, ρ)
    σ = 10.0
    β = 8 / 3
    u0 = [x₀, 1.0, 0.0]
    tspan = (0.0, 100.0)

    function lorenz!(du, u, p, t)
        x, y, z = u
        du[1] = σ * (y - x)
        du[2] = x * (ρ - z) - y
        du[3] = x * y - β * z
        return nothing
    end

    prob = ODEProblem(lorenz!, u0, tspan)
    return solve(prob)
end

function lorenz_plot()
    sol14_1 = lorenz_solver(0.0, 14.0)
    sol14_2 = lorenz_solver(1.0e-5, 14.0)

    sol28_1 = lorenz_solver(0.0, 28.0)
    sol28_2 = lorenz_solver(1.0e-5, 28.0)

    # Common grid for top four plots
    t = 0.0:0.1:100.0

    u14_1 = sol14_1(t).u
    u14_2 = sol14_2(t).u
    u28_1 = sol28_1(t).u
    u28_2 = sol28_2(t).u

    x14_1 = getindex.(u14_1, 1)
    x14_2 = getindex.(u14_2, 1)
    x28_1 = getindex.(u28_1, 1)
    x28_2 = getindex.(u28_2, 1)

    diff14 = abs.(x14_1 .- x14_2)
    diff28 = abs.(x28_1 .- x28_2)

    # Dense grid for smoother phase plots
    tphase = 0.0:0.01:100.0
    uph14 = sol14_1(tphase).u
    uph28 = sol28_1(tphase).u

    x14_phase = getindex.(uph14, 1)
    z14_phase = getindex.(uph14, 3)

    x28_phase = getindex.(uph28, 1)
    z28_phase = getindex.(uph28, 3)

    # Top-left: show y labels, hide x numbers
    p1 = plot(
        t, x14_1;
        title = L"\rho = 14.0",
        ylabel = L"x_1,\ x_2",
        legend = false,
        xformatter = blankfmt,
        yformatter = tickfmt
    )
    plot!(p1, t, x14_2)

    # Top-right: hide x numbers and y numbers, keep tick marks
    p2 = plot(
        t, x28_1;
        title = L"\rho = 28.0",
        legend = false,
        xformatter = blankfmt,
        yformatter = blankfmt
    )
    plot!(p2, t, x28_2)

    # Middle-left: show both axes
    p3 = plot(
        t, diff14;
        xlabel = L"t",
        ylabel = L"|x_1 - x_2|",
        legend = false,
        xformatter = tickfmt,
        yformatter = tickfmt
    )

    # Middle-right: show x numbers, hide y numbers
    p4 = plot(
        t, diff28;
        xlabel = L"t",
        legend = false,
        xformatter = tickfmt,
        yformatter = blankfmt
    )

    # Bottom-left: show both axes
    p5 = plot(
        x14_phase, z14_phase;
        xlabel = L"x_1",
        ylabel = L"z_1",
        legend = false,
        xformatter = tickfmt,
        yformatter = tickfmt
    )

    # Bottom-right: show x numbers, hide y numbers
    p6 = plot(
        x28_phase, z28_phase;
        xlabel = L"x_1",
        legend = false,
        xformatter = tickfmt,
        yformatter = blankfmt
    )

    topblock = plot(p1, p2, p3, p4; layout = (2, 2), link = :both)
    bottomblock = plot(p5, p6; layout = (1, 2), link = :both)

    return plot(
        topblock, bottomblock;
        layout = grid(2, 1, heights = [2/3, 1/3])
    )
end

export lorenz_solver, lorenz_plot

end
