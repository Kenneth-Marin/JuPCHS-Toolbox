"""
# Operations that preserves structure
# function Parallel(S1::typePCHS, S2::typePCHS)
# function Feedbacks(S1::typePCHS, S2::typePCHS)
"""

function check_repeated_label(label::Vector)
    n = length(label)
    for k = 1:n
        xl = label[k]
        r = findall(label.==xl)
        s = length(r)
        if s > 1 # hay variables repetidas
           for m = 2:s
               label[r[m]] = label[r[m]]*"p" 
           end 
        end
    end
    return label
end

function Parallel(S1::typePCHS, S2::typePCHS)
    println("Parallel connection")
    S3 = deepcopy(S1)
    if S1.nu != S2.nu
        print_warning("The inputs are not consistents") 
     end 
    if S1.dt != S2.dt
        S3.dt = minimum([S1.dt S2.dt])  
         print_warning("Different discretization time. Using dt = "*string(S3.dt))
    end
    S3.nx = S1.nx + S2.nx 
    p1 = 1:S1.nx
    p2 = (S1.nx+1):S3.nx
    S3.name = "["*S1.name*"+"*S2.name*"]"
    S3.xlabel = check_repeated_label([S1.xlabel;S2.xlabel])
    for k = 1:S1.nu
        if S1.ulabel[k]!=S2.ulabel[k]
           S3.ulabel[k] = S1.ulabel[k] * S2.ulabel[k]
        end
        if S1.ylabel[k]!=S2.ylabel[k]
           S3.ylabel[k] = S1.ylabel[k] * S2.ylabel[k]
        end
    end
    S3.H = x -> S1.H(x[p1]) + S2.H(x[p2])
    S3.dH = x -> [S1.dH(x[p1]); S2.dH(x[p2])]
    S3.J = x -> [S1.J(x[p1],u[p1]) zeros(S1.nx,S2.nx); zeros(S2.nx,S1.nx) S2.J(x[p2],u[p2])]
    S3.R = x -> [S1.R(x[p1]) zeros(S1.nx,S2.nx); zeros(S2.nx,S1.nx) S2.R(x[p2])]
    S3.G = x -> [S1.G(x[p1]); S2.G(x[p1])]
    return S3
end

function RepeatParallel(m::Int, S1::typePCHS)
    S3 = deepcopy(S1)
    if m > 1
        for k = 2:m
            S3 = S3 + S1
        end
    end
    return S3
end

function Feedback(S1::typePCHS, S2::typePCHS)
    S3 = deepcopy(S1)
    println("Feedback connection")
    if S1.nu != S2.nu
       print_warning("The inputs are not consistents") 
     end 
    if S1.dt != S2.dt
       S3.dt = minimum([S1.dt S2.dt])  
       print_warning("Different discretization time. Using dt = "*string(S3.dt))
    end
    S3.nx = S1.nx + S2.nx
    S3.nu = S1.nu + S2.nu
    S3.name = "["*S1.name*"*"*S2.name*"]"
    S3.xlabel = check_repeated_label([S1.xlabel;S2.xlabel])
    S3.ulabel = check_repeated_label([S1.ulabel;S2.ulabel])
    S3.ylabel = check_repeated_label([S1.ylabel;S2.ylabel])
    p1 = 1:S1.nx
    p2 = (S1.nx+1):S3.nx
    S3.H = x -> S1.H(x[p1]) + S2.H(x[p2]) 
    S3.dH = x -> [S1.dH(x[p1]); S2.dH(x[p2])]
    S3.R = x -> [S1.R(x[p1]) zeros(S1.nx,S2.nx); zeros(S2.nx,S1.nx) S2.R(x[p2])]
    S3.J = x -> [S1.J(x[p1],u[p2]) -S1.G(x[p1])*transpose(S2.G(x[p2])); S2.G(x[p2])*transpose(S1.G(x[p1])) S2.J(x[p2],u[p2])]
    S3.G = x -> [S1.G(x[p1]) zeros(S1.nx,S1.nu); zeros(S2.nx,S2.nu) S2.G(x[p2])]
    return S3
end

function MatrixTransformation(S1::typePCHS, M::Matrix)
    """
    S1: typePCHS
    M: Matrix to transform the system
    """
    S3 = deepcopy(S1)
    println("Matrix transformation")
    if S1.nu != size(M,1)
        print_error("The inputs are not consistents with the matrix size") 
        return nothing
    end
    S3.nx = S1.nx
    S3.nu = S1.nu
    S3.name = "[M*"*S1.name*"M^T]"
    # Realizar la multiplicación
    S3.H = x -> S1.H(x)
    S3.dH = x -> S1.dH(x)
    S3.R = x -> S1.R(x)
    S3.J = x -> S1.J(x,u)
    S3.G = x -> S1.G(x)*M

    return S3
end

# Incremental model (Aun no funciona )
function IncrementalModel(;PCHS=nothing,xp=nothing, up=nothing)
    if PCHS===nothing
        print_error("PCHS is empty")
        return nothing
    end
    if xp===nothing
        print_error("xp is empty")
        return nothing
    end
    if up===nothing
        print_error("up is empty")
        return nothing
    end    
    println("Incremental model")
    println("\t∆x' = (J(∆x)-R(∆x))∇H(∆x) + G(∆x)u(∆x)")

    PCHS.H = x -> PCHS.H(x-xp)
    PCHS.dH = x -> PCHS.dH(x-xp)
    PCHS.R = x -> PCHS.R(x-xp)
    PCHS.J = x -> PCHS.J(x-xp,u)
    PCHS.G = x -> PCHS.G(x-xp)
    #PCHS.u = u -> u - up
end 

# overload + operation
import Base.+
(+)(S1::typePCHS, S2::typePCHS) = Parallel(S1,S2)

# overload * operation
import Base.*
(*)(S1::typePCHS, S2::typePCHS) = Feedback(S1,S2)
(*)(m::Int, S1::typePCHS) = RepeatParallel(m,S1)
(*)(S1::typePCHS, m::Int) = RepeatParallel(m,S1)
(*)(M::Matrix, S1::typePCHS) = MatrixTransformation(S1,M)
(*)(S1::typePCHS, M::Matrix) = MatrixTransformation(S1,M)
