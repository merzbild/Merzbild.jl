@testset "1-D Poisson solver" begin
    ε0 = Merzbild.eps_0

    node_coords(grid) = [(j - 1) * grid.Δx for j in 1:grid.n_cells+1]

    L = 0.05
    nx = 40
    grid = Grid1DUniform(L, nx)
    x = node_coords(grid)

    # ρ = 0, Dirichlet-Dirichlet: linear potential, constant field in every node
    V = 12.0
    poisson_solver = PoissonSolver1DUniform(grid, DirichletFieldBC1D(0.0), DirichletFieldBC1D(V))
    field_props = ElectrostaticFieldProps(grid)

    @test field_props.n_nodes == nx + 1
    @test poisson_solver.n_unknowns == nx - 1
    @test poisson_solver.node_offset == 1

    solve_poisson!(poisson_solver, field_props)

    @test maximum(abs.(field_props.potential .- V .* x ./ L)) < 1e-13 * V
    @test maximum(abs.(field_props.electric_field .+ V / L)) < 1e-11 * V / L
    @test field_props.net_charge_density == 0.0

    # a field props instance created from the number of nodes is interchangeable with one created from the grid
    field_props_n = ElectrostaticFieldProps(nx + 1)
    solve_poisson!(poisson_solver, field_props_n)
    @test field_props_n.potential == field_props.potential
    @test field_props_n.electric_field == field_props.electric_field

    # mutating the BC value between solves changes the solution accordingly
    poisson_solver.bc_right.ϕ = -V
    solve_poisson!(poisson_solver, field_props)
    @test maximum(abs.(field_props.potential .+ V .* x ./ L)) < 1e-13 * V
    poisson_solver.bc_right.ϕ = V

    # ρ = const, Dirichlet-Dirichlet: the parabolic solution is exact for a 2nd-order stencil
    ρ0 = 3e-7
    poisson_solver = PoissonSolver1DUniform(grid, DirichletFieldBC1D(0.0), DirichletFieldBC1D(0.0))
    field_props = ElectrostaticFieldProps(grid)
    fill!(field_props.charge_density, ρ0)
    solve_poisson!(poisson_solver, field_props)

    ϕ_exact = ρ0 .* x .* (L .- x) ./ (2 * ε0)
    @test maximum(abs.(field_props.potential .- ϕ_exact)) < 1e-13 * maximum(ϕ_exact)

    # ρ = 0, Neumann-Dirichlet: ϕ = E_L (L - x), E ≡ E_L
    E_L = 500.0
    poisson_solver = PoissonSolver1DUniform(grid, NeumannFieldBC1D(E_L), DirichletFieldBC1D(0.0))
    field_props = ElectrostaticFieldProps(grid)

    @test poisson_solver.n_unknowns == nx
    @test poisson_solver.node_offset == 0

    solve_poisson!(poisson_solver, field_props)

    @test maximum(abs.(field_props.potential .- E_L .* (L .- x))) < 1e-13 * E_L * L
    @test maximum(abs.(field_props.electric_field .- E_L)) < 1e-11 * E_L
    @test field_props.electric_field[1] == E_L

    # ρ = 0, Dirichlet-Neumann: ϕ = ϕ_L - E_R x
    E_R = -320.0
    ϕ_L = 7.0
    poisson_solver = PoissonSolver1DUniform(grid, DirichletFieldBC1D(ϕ_L), NeumannFieldBC1D(E_R))
    field_props = ElectrostaticFieldProps(grid)

    @test poisson_solver.n_unknowns == nx
    @test poisson_solver.node_offset == 1

    solve_poisson!(poisson_solver, field_props)

    @test maximum(abs.(field_props.potential .- (ϕ_L .- E_R .* x))) < 1e-12 * abs(E_R) * L
    @test maximum(abs.(field_props.electric_field .- E_R)) < 1e-11 * abs(E_R)
    @test field_props.electric_field[nx+1] == E_R

    # ρ = const ≠ 0, Neumann-Dirichlet: checks the halved boundary row against the ×2 density normalization
    ϕ_R = -2.0
    poisson_solver = PoissonSolver1DUniform(grid, NeumannFieldBC1D(E_L), DirichletFieldBC1D(ϕ_R))
    field_props = ElectrostaticFieldProps(grid)
    fill!(field_props.charge_density, ρ0)
    solve_poisson!(poisson_solver, field_props)

    C = ϕ_R + ρ0 * L^2 / (2 * ε0) + E_L * L
    ϕ_exact = -ρ0 .* x.^2 ./ (2 * ε0) .- E_L .* x .+ C
    @test maximum(abs.(field_props.potential .- ϕ_exact)) < 1e-13 * maximum(abs.(ϕ_exact))

    # periodic with ρ = const: the mean subtraction leaves a vanishing RHS
    poisson_solver = PoissonSolver1DUniform(grid, PeriodicFieldBC1D(), PeriodicFieldBC1D())
    field_props = ElectrostaticFieldProps(grid)
    fill!(field_props.charge_density, ρ0)

    @test poisson_solver.n_unknowns == nx - 1
    @test poisson_solver.node_offset == 0

    solve_poisson!(poisson_solver, field_props)

    @test maximum(abs.(field_props.potential)) < 1e-13 * ρ0 * L^2 / ε0
    @test maximum(abs.(field_props.electric_field)) < 1e-13 * ρ0 * L / ε0
    @test abs(field_props.net_charge_density - ρ0) < 1e-14 * ρ0

    # convergence: Dirichlet-Dirichlet, ρ = ρ0 sin(π x / L)
    errors = Float64[]
    for nxc in [20, 40, 80]
        g = Grid1DUniform(L, nxc)
        xx = node_coords(g)
        ps = PoissonSolver1DUniform(g, DirichletFieldBC1D(0.0), DirichletFieldBC1D(0.0))
        fp = ElectrostaticFieldProps(g)

        for j in 1:nxc+1
            fp.charge_density[j] = ρ0 * sin(π * xx[j] / L)
        end
        solve_poisson!(ps, fp)

        exact = ρ0 * L^2 / (π^2 * ε0) .* sin.(π .* xx ./ L)
        push!(errors, maximum(abs.(fp.potential .- exact)))
    end
    @test abs(errors[1] / errors[2] - 4.0) < 0.1
    @test abs(errors[2] / errors[3] - 4.0) < 0.1

    # convergence: periodic, ρ = ρ0 sin(2π k x / L); the exact solution has zero mean, so the gauge is checked as well
    k_wave = 2
    errors = Float64[]
    for nxc in [40, 80, 160]
        g = Grid1DUniform(L, nxc)
        xx = node_coords(g)
        ps = PoissonSolver1DUniform(g, PeriodicFieldBC1D(), PeriodicFieldBC1D())
        fp = ElectrostaticFieldProps(g)

        for j in 1:nxc+1
            fp.charge_density[j] = ρ0 * sin(2π * k_wave * xx[j] / L)
        end
        solve_poisson!(ps, fp)

        exact = ρ0 * L^2 / (4 * π^2 * k_wave^2 * ε0) .* sin.(2π * k_wave .* xx ./ L)
        push!(errors, maximum(abs.(fp.potential .- exact)))

        @test abs(sum(fp.potential[1:nxc])) < 1e-12 * maximum(abs.(exact))
        @test fp.potential[nxc+1] == fp.potential[1]
        @test fp.electric_field[nxc+1] == fp.electric_field[1]
    end
    @test abs(errors[1] / errors[2] - 4.0) < 0.1
    @test abs(errors[2] / errors[3] - 4.0) < 0.1

    # the Dirichlet-Dirichlet Poisson FD matrix factorization is analytic: cp[i] = -i/(i+1), inv_m[i] = Δx² i/(i+1)
    g = Grid1DUniform(L, 30)
    ps = PoissonSolver1DUniform(g, DirichletFieldBC1D(0.0), DirichletFieldBC1D(0.0))
    n_unknowns = ps.n_unknowns

    for i in 1:n_unknowns-1
        @test abs(ps.cp[i] + i / (i + 1)) < 4 * eps()
    end
    @test ps.cp[n_unknowns] == 0.0  # c[n_unknowns] lies outside of the matrix and is stored as 0
    for i in 1:n_unknowns
        @test abs(ps.inv_m[i] / (g.Δx^2 * i / (i + 1)) - 1.0) < 4 * eps()
    end

    # residual check A ϕ = rhs in every admissible BC combination, using the stored diagonals
    for (bc_left, bc_right) in [(DirichletFieldBC1D(1.0), DirichletFieldBC1D(-3.0)),
                                (NeumannFieldBC1D(120.0), DirichletFieldBC1D(2.0)),
                                (DirichletFieldBC1D(2.0), NeumannFieldBC1D(-40.0)),
                                (PeriodicFieldBC1D(), PeriodicFieldBC1D())]
        g = Grid1DUniform(L, 25)
        ps = PoissonSolver1DUniform(g, bc_left, bc_right)
        fp = ElectrostaticFieldProps(g)

        for j in 1:fp.n_nodes
            fp.charge_density[j] = ρ0 * sin(3.0 * j)
        end
        solve_poisson!(ps, fp)

        n = ps.n_unknowns
        residual = 0.0
        for i in 1:n
            row = ps.b[i] * ps.x[i]
            if i > 1
                row += ps.a[i] * ps.x[i-1]
            end
            if i < n
                row += ps.c[i] * ps.x[i+1]
            end
            residual = max(residual, abs(row - ps.rhs[i]))
        end
        @test residual < 1e-12 * maximum(abs.(ps.rhs))

        # the values of the potential in the unknown nodes are the solution of the tridiagonal system
        # (up to the constant gauge shift applied in the periodic case)
        if !isa(bc_left, PeriodicFieldBC1D)
            for i in 1:n
                @test fp.potential[ps.node_offset + i] == ps.x[i]
            end
        end
    end

    # periodic: the residual of the dropped row is ~roundoff after the mean subtraction
    g = Grid1DUniform(L, 25)
    ps = PoissonSolver1DUniform(g, PeriodicFieldBC1D(), PeriodicFieldBC1D())
    fp = ElectrostaticFieldProps(g)
    for j in 1:fp.n_nodes
        fp.charge_density[j] = ρ0 * sin(3.0 * j)
    end
    fp.charge_density[fp.n_nodes] = fp.charge_density[1]
    solve_poisson!(ps, fp)

    n_cells = g.n_cells
    dropped_row = (-fp.potential[n_cells-1] + 2 * fp.potential[n_cells] - fp.potential[1]) * g.inv_Δx^2
    dropped_rhs = (fp.charge_density[n_cells] - fp.net_charge_density) / ε0
    @test abs(dropped_row - dropped_rhs) < 1e-12 * maximum(abs.(ps.rhs))
    @test abs(sum(fp.potential[1:n_cells])) < 1e-12 * maximum(abs.(fp.potential))
    @test fp.potential[fp.n_nodes] == fp.potential[1]
    @test fp.electric_field[fp.n_nodes] == fp.electric_field[1]

    # the pinned node is written back at every solve and not left over from the previous one
    ϕ_pinned_first = fp.potential[n_cells]
    fp.potential[n_cells] = 1e5
    solve_poisson!(ps, fp)
    @test abs(fp.potential[n_cells] - ϕ_pinned_first) < 1e-12 * maximum(abs.(fp.potential))

    # two consecutive solves with different RHS give the same result as two freshly constructed solvers
    g = Grid1DUniform(L, 32)
    ps = PoissonSolver1DUniform(g, DirichletFieldBC1D(3.0), NeumannFieldBC1D(-11.0))
    fp = ElectrostaticFieldProps(g)

    ps_fresh1 = PoissonSolver1DUniform(g, DirichletFieldBC1D(3.0), NeumannFieldBC1D(-11.0))
    fp_fresh1 = ElectrostaticFieldProps(g)
    ps_fresh2 = PoissonSolver1DUniform(g, DirichletFieldBC1D(3.0), NeumannFieldBC1D(-11.0))
    fp_fresh2 = ElectrostaticFieldProps(g)

    for j in 1:fp.n_nodes
        fp.charge_density[j] = ρ0 * cos(0.7 * j)
        fp_fresh1.charge_density[j] = fp.charge_density[j]
    end
    solve_poisson!(ps, fp)
    solve_poisson!(ps_fresh1, fp_fresh1)
    @test fp.potential == fp_fresh1.potential

    for j in 1:fp.n_nodes
        fp.charge_density[j] = ρ0 * sin(1.3 * j)
        fp_fresh2.charge_density[j] = fp.charge_density[j]
    end
    solve_poisson!(ps, fp)
    solve_poisson!(ps_fresh2, fp_fresh2)
    @test fp.potential == fp_fresh2.potential
    @test fp.electric_field == fp_fresh2.electric_field

    # solving twice on the same charge density changes nothing
    ϕ_first = copy(fp.potential)
    solve_poisson!(ps, fp)
    @test fp.potential == ϕ_first

    # inadmissible BC combinations
    @test_throws ErrorException PoissonSolver1DUniform(grid, NeumannFieldBC1D(0.0), NeumannFieldBC1D(0.0))
    @test_throws ErrorException PoissonSolver1DUniform(grid, PeriodicFieldBC1D(), DirichletFieldBC1D(0.0))
    @test_throws ErrorException PoissonSolver1DUniform(grid, DirichletFieldBC1D(0.0), PeriodicFieldBC1D())
    @test_throws ErrorException PoissonSolver1DUniform(grid, PeriodicFieldBC1D(), NeumannFieldBC1D(0.0))

    # a grid too coarse for the chosen BCs
    @test_throws ErrorException PoissonSolver1DUniform(Grid1DUniform(L, 1),
                                                       DirichletFieldBC1D(0.0), DirichletFieldBC1D(0.0))

    # clearing the field properties
    fp_clear = ElectrostaticFieldProps(grid)
    fill!(fp_clear.charge_density, 1.0)
    fill!(fp_clear.potential, 2.0)
    fill!(fp_clear.electric_field, 3.0)
    fp_clear.net_charge_density = 4.0

    clear_charge_density!(fp_clear)
    @test maximum(abs.(fp_clear.charge_density)) == 0.0
    @test fp_clear.net_charge_density == 0.0
    @test minimum(fp_clear.potential) == 2.0
    @test minimum(fp_clear.electric_field) == 3.0

    fill!(fp_clear.charge_density, 1.0)
    clear_props!(fp_clear)
    @test maximum(abs.(fp_clear.charge_density)) == 0.0
    @test maximum(abs.(fp_clear.potential)) == 0.0
    @test maximum(abs.(fp_clear.electric_field)) == 0.0
    @test fp_clear.net_charge_density == 0.0
end
