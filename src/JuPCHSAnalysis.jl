"""
Tools for analysis
# function Heaviside(t)
# function AdmissibleEquilibrium(;PCHS=nothing,u=nothing,x=nothing)
# function Hmin(;PCHS=nothing,xini=nothing)
"""


# Step variable for the optimization methods
function line_search(f, g, x, p; α_init=1.0, c1=1e-4, tau=0.5)
    α = α_init
    while f(x + α * p) > f(x) + c1 * α * (g(x)' * p)
        α *= tau
    end
    return α
end

function convergence_plot(x; xlabel="Iterations", ylabel="Norm of gradient", title="", label="", add=false)
    if add
        pl = plot!(x, title=title, xlabel=xlabel, ylabel=ylabel, label=label)
    else
        pl = plot(x, title=title, xlabel=xlabel, ylabel=ylabel, label=label)
    end 
    return pl
end

# No-nlinear conjugate gradient method  
function nonlinear_conjugate_gradient(f, g, x0; tol=EPSILON, max_iter=MAXITER)
    x = x0
    grad = g(x)
    p = -grad
    grad_norm = zeros(max_iter,1)
    for k = 1:max_iter
        grad_norm[k] = norm(grad)
        if grad_norm[k] < tol
            break
        end
        α = line_search(f, g, x, p)
        x_new = x + α * p
        grad_new = g(x_new)
        β = max(0, (grad_new' * (grad_new - grad)) / (grad' * grad))
        p = -grad_new + β * p
        x = x_new
        grad = grad_new
    end
    return x, grad_norm
end

# Gradient method 
function gradient(f, g, x0; tol=EPSILON, max_iter=MAXITER, fixed=nothing)
    """
    Gradient method
    f: function to minimize
    g: gradient of f
    x0: initial point
    tol: tolerance
    max_iter: maximum number of iterations
    fixed: fixed step size
    """
    if fixed === nothing
        bool = true
    elseif fixed isa Number
        bool = false
    else 
        error("Invalid fixed step size. It will be a 0.1")
        fixed = 0.1
        bool = false
    end

    x = x0                         # Initial point
    grad = g(x)                    # Gradient at x
    grad_norm = zeros(max_iter,1)  # Norm of the gradient
    
    for k = 1:max_iter
        if bool
            α = line_search(f, g, x, -grad)
        else
            α = fixed
        end
        x = x - α * grad
        grad = g(x)
        grad_norm[k] = norm(grad)
        if grad_norm[k] < tol
            break
        end
    end
    return x, grad_norm
end

# Newton method
function newton(f, df, x0; tol=1e-6, max_iter=100)
    """
    f: function to minimize
    df: gradient of f
    x0: initial point
    tol: tolerance
    max_iter: maximum number of iterations
    """
    hes(x) = ForwardDiff.hessian(f, x)   # Hessian of f

    # Initial point
    x = x0
    grad = df(x)
    hess = hes(x)
    grad_norm = zeros(max_iter,1)

    for k = 1:max_iter
        grad_norm[k] = norm(grad)
        if grad_norm[k] < tol
            break
        end
        x = x - hess \ grad
        grad = df(x)    # Gradient evaluated
        hess = hes(x)   # Hessian evaluated
    end
    return x, grad_norm
end

# Elipsoid method
function elipsoid(f, df, x0; tol=1e-6, max_iter=100)
    """
    f: function to minimize
    df: gradient of f
    x0: initial point
    tol: tolerance
    max_iter: maximum number of iterations
    """
    n = length(x0)
    B = Matrix(I, n, n)  # Init B 
    c = x0               # Init center
    grad_norm = zeros(max_iter,1)

    for k in 1:max_iter
        α = df(c)         # Gradient evaluated
        grad_norm[k] = norm(α)
        if grad_norm[k] < tol
            break
        end   
        # vector p
        αTBBTα = sqrt(α' * B * B' * α)
        p = (1 / αTBBTα) * B' * α
        # Uploading center c+
        c = c - (1 / (n + 1)) * B * p
        # Uploading B+
        B = (n / (sqrt(n^2 - 1))) * B + ((n / (n + 1)) - (n / sqrt(n^2 - 1))) * (B * p) * p'
    end
    
    return c, grad_norm
end

# Fixed point method for equilibrium point (f(x,u)=0)
function fixed_point_(PCHS, x0, u; tol=1e-6, max_iter=100)    
    # Initial conditions
    x = x0
    ut = u(x,0) 
    F(x,ur) = 0.1*(PCHS.J(x,u)-PCHS.R(x))*PCHS.dH(x)+PCHS.G(x)*ur
    # Metodo de punto fijo
    for k = 1:max_iter
        xn = x + F(x,ut)
        err = norm(xn-x)
        if err < tol
            break
        end
        # Actualizar
        x = xn
        ut = u(x,(k-1))
    end 
    return x, ut
end

# CASE i: Given a desired x, it will optimize
function kkt_conditions_x(PCHS, xd, v::AbstractVector)
    nx = PCHS.nx
    nu = PCHS.nu
    # Only one argument
    x = v[1:nx]
    u = v[nx+1:nx+nu]
    l = v[nx+nu+1:2*nx+nu]
    
    # Hamiltonian structure
    F(x, u) = (PCHS.J(x,u) - PCHS.R(x)) * PCHS.dH(x) + PCHS.G(x) * u
    
    # Partial derivatives
    dfdx = ForwardDiff.jacobian(x -> F(x, u), x)
    dfdu = ForwardDiff.jacobian(u -> F(x, u), u)
    
    # KKT conditions:
    eq1 = x - xd + dfdx' * l
    eq2 = dfdu' * l
    eq3 = F(x, u)
    
    # Only one equation
    return [eq1; eq2; eq3]
end

# CASE ii: Given a desired u, it will optimize
function kkt_conditions_u(PCHS, ud, v::AbstractVector)
    nx = PCHS.nx
    nu = PCHS.nu
    # Only one argument
    x = v[1:nx]
    u = v[nx+1:nx+nu]
    l = v[nx+nu+1:2*nx+nu]
    
    # Hamiltonian structure
    F(x, u) = (PCHS.J(x,u) - PCHS.R(x)) * PCHS.dH(x) + PCHS.G(x) * u
    
    # Partial derivatives
    dfdx = ForwardDiff.jacobian(x -> F(x, u), x)
    dfdu = ForwardDiff.jacobian(u -> F(x, u), u)
    
    # KKT conditions:
    eq1 = dfdx' * l
    eq2 = u - ud + dfdu' * l
    eq3 = F(x, u)
    
    # Only one equation
    return [eq1; eq2; eq3]
end

# Newton method for equilibrium point (f(x,u)=0)
function newton_(f, x0; tol=1e-6, max_iter=100)
    x = x0
    for iter in 1:max_iter
        fx = f(x)                           # Function evaluated
        J = ForwardDiff.jacobian(f, x)      # Jacobian evaluated
        delta = J \ fx                      # Newton step
        x_new = x - delta                   
        x = x_new                           # Update
        if norm(x_new - x) < tol
            break
        end
    end
    return x, max_iter
end

# _________________________________________________________________________
#                           External functions
# _________________________________________________________________________

function Heaviside(t)
    return 0.5 * (sign(t) + 1)
end

function AdmissibleEquilibrium(;PCHS=nothing,u=nothing,x=nothing)
    if PCHS===nothing
        print_warning("Empty PCHS")
        PCHS =  build_PHS()
    end
    if u===nothing
        print_warning("Assuming input u=0")
        u = zeros(PCHS.nu)
    end
    if x===nothing
        print_warning("Assuming states x=0")
        x = zeros(PCHS.nx)
    end
    A = (PCHS.J(x,u)-PCHS.R(x))
    w = (A*PCHS.dH(x) + PCHS.G(x)*u)/norm(A)
    r = norm(w)
    rs = true
    if r > EPSILON
        print_warning("It is not an admissible equilibrium")
        println("|Δ| =",r, ", dx=", w)
        rs = false
    end
    return rs
end

# Find the minimum of the Hamiltonian
function Hmin(;PCHS=nothing,xini=nothing, method=nothing, all=false, alpha=nothing, tol=EPSILON, max_iter=MAXITER)
    """
    PCHS : port-controlled Hamiltonian system
    xini : initial point
    method : optimization method
    all : all methods
    alpha : fixed step size
    """
    if PCHS === nothing
        print_warning("Empty PCHS")
        PCHS =  build_PHS()
    end
    if xini === nothing
        xini = zeros(PCHS.nx)
    end
    if method === nothing && all == false
        print_warning("Not selected a method. It will be the newton method")
        method = "newton"
    end
    if all
        println("Optimizing by all methods")
    else
        println("Optimizing by the ", method, " method")
    end

    plot()   # Reset the graph
    if method in ["nonlinear gradient", "nonlinear conjugate gradient"] || all
        x_opt, grad_norm = nonlinear_conjugate_gradient(PCHS.H, PCHS.dH, xini, tol=tol, max_iter=max_iter)
        pl = convergence_plot(grad_norm, title="Convergence plot", label="nonlinear conjugate gradient", add=all)
    end 
    if method in ["gradient", "Gradient"] || all
        x_opt, grad_norm = gradient(PCHS.H, PCHS.dH, xini, tol=tol, max_iter=max_iter)
        pl = convergence_plot(grad_norm, title="Convergence plot", label="gradient", add=all)
    end
    if method in ["fixed gradient", "Fixed Gradient"] || all
        if alpha === nothing 
            if all == false
                print_warning("Not selected a fixed step size. It will be 0.1")
            end 
            alpha = 0.1
        end
        label = "fixed gradient (t = $(alpha))"
        x_opt, grad_norm = gradient(PCHS.H, PCHS.dH, xini, fixed=alpha, tol=tol, max_iter=max_iter)
        pl = convergence_plot(grad_norm, title="Convergence plot", label=label, add=all)
    end
    if method in ["elipsoid", "Elipsoid"] || all
        x_opt, grad_norm = elipsoid(PCHS.H, PCHS.dH, xini, tol=tol, max_iter=max_iter)
        pl = convergence_plot(grad_norm, title="Convergence plot", label="Elipsoid", add=all)
    end
    if method in ["newton", "Newton"] || all
        x_opt, grad_norm = newton(PCHS.H, PCHS.dH, xini, tol=tol, max_iter=max_iter)
        pl = convergence_plot(grad_norm, title="Convergence plot", label="Newton", add=all)
    end

    # Convergence plot
    display(pl)
    hm = PCHS.H(x_opt)
    println("Hamiltonian minimum value =", hm)

    return x_opt, hm
end

# Find the equilibrium point of the system (f(x,u)=0)
function equilibrium_point(;PCHS=nothing,xini=nothing, uini=nothing,
    xd=nothing,ud=nothing, method=nothing, tol=EPSILON, max_iter=MAXITER)
    """
    PCHS : port-controlled Hamiltonian system
    xini : initial point
    uini : initial input
    xd : desired x
    ud : desired u
    method : optimization method
    """
    if PCHS === nothing
        print_warning("Empty PCHS")
        PCHS =  build_PHS()
    end
    if xini === nothing
        xini = zeros(PCHS.nx)
    end
    if uini === nothing
        uini = zeros(PCHS.nu)
    end
    if xd === nothing && ud === nothing
        print_warning("Empty desired x and u. X desired will be assumed as the initial point")
        xd = xini
    end
    if method === nothing
        print_warning("Not selected a method. It will be the newton method")
        method = "newton"
    end

    if method in ["fixed point", "Fixed Point"]
        x_opt, u_opt = fixed_point_(PCHS, xini, uini, tol=tol, max_iter=max_iter)
    elseif method in ["Newton", "newton"]
        if xd === nothing
            # Conditions kkt for u
            f33 = v -> kkt_conditions_u(PCHS, ud, v)
            # Only one argument v = [x; u; l]
            v0 = [xini; uini; zeros(PCHS.nx)]
            # Optimization
            vopt, iter = newton_(f33, v0, tol=tol, max_iter=max_iter)
            x_opt = vopt[1:PCHS.nx]
            u_opt = vopt[PCHS.nx+1:PCHS.nx+PCHS.nu]
        else
            # Conditions kkt for x
            f33 = v -> kkt_conditions_x(PCHS, xd, v)
            # Only one argument v = [x; u; l]
            v0 = [xini; uini; zeros(PCHS.nx)]
            # Optimization
            vopt, iter = newton_(f33, v0, tol=tol, max_iter=max_iter)
            x_opt = vopt[1:PCHS.nx]
            u_opt = vopt[PCHS.nx+1:PCHS.nx+PCHS.nu]
        end
    end
    return x_opt, u_opt
end

# Numerical verification of Lipschitz continuity
function islipschitz(;f=nothing,inf=nothing,sup=nothing, nx=nothing, nt=nothing, tol=0.01)
    """
    f : function to verify
    inf : lower bound
    sup : upper bound
    nx : dimension
    nt : number of iterations
    """
    if f === nothing || nx === nothing
        print_error("Empty function or dimension")
        return nothing
    end
    if nt === nothing
        nt = MAXITER
    end
    # Calcula el 70% de las iteraciones 
    n_rev = round(Int, nt*0.7)
    # Inicializar L
    L = 0
    iPlot = 1
    # Iterar sobre un conjunto de puntos aleatorios
    iter = zeros(nt)
    for i in 1:nt
        # Comprobar que es Lipschitz
        if inf === nothing && sup === nothing
            x = rand(nx)
            y = rand(nx)
        else
            if inf === nothing 
                inf = -100
            end
            if sup === nothing 
                sup = 100
            end
            x = (inf + (sup-inf)).*rand(nx)
            y = (inf + (sup-inf)).*rand(nx)
        end        
        # Aproximación de la constante de Lipschitz
        Lx = norm(f(x) - f(y)) / norm(x - y)
        # Actualizar L
        L = max(Lx,L)
        iter[i] = L
        
        # Revisar que la constante permanezca constante
        if i >= n_rev
            if iter[i-1] < L - tol || iter[i-1] > L + tol
                print_warning("It could not find a constant Lipschitz")
                L = nothing
                iPlot = 0
                break
            end
        end 
    end

    if iPlot == 1
        # Graficar
        plt = plot(iter, label="Lipschitz")
        plt = plot!(xlabel="Iteraciones", ylabel="Lipschitz")
        plt = plot!(title="Convergencia de la constante de Lipschitz")
        display(plt)
    end

    return L
end 
