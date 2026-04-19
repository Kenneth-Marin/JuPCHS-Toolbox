"""
# Display and plot functions #
# function Base.show(io::IO, p::typePCHS) #
# function print_warning(s::String) #
# function plot_results(;PCHS=nothing,data=nothing,title="PCHS system") #        
"""

# Printing the object in the terminal
function Base.show(io::IO, p::typePCHS)
    printstyled(io,"📌 $(p.name) (Port Controlled Hamiltonian System)\n", color=:blue)
    println("\tx' = (J(x)-R(x))∇H(x) + G(x)u(x)")    
    if (p.nx<4)&&(p.nx>0)
        s = "State variables:\t["
        for k = 1:p.nx-1
            s = s*p.xlabel[k]*", "
        end
        s = s*p.xlabel[p.nx]*"]"
        println(io,s)
    else
        println(io,"State variables:\t[",p.xlabel[1],"...",p.xlabel[p.nx],"]")
    end    
    if (p.nu<4)&&(p.nu>0)
        s = "Inputs:         \t["
        for k = 1:p.nu-1
            s = s*p.ulabel[k]*", "
        end
        s = s*p.ulabel[p.nu]*"]"
        println(io,s)
    else
        println(io,"Inputs:         \t[",p.ulabel[1],"...",p.ulabel[p.nu],"]\t")
    end    
    if (p.nu<4)&&(p.nu>0)
        s = "Outputs:         \t["
        for k = 1:p.nu-1
            s = s*p.ylabel[k]*", "
        end
        s = s*p.ylabel[p.nu]*"]"
        println(io,s)
    else
        println(io,"Outputs:         \t[",p.ylabel[1],"...",p.ylabel[p.nu],"]\t")
    end    
    println(io,"Discretization time:\t$(p.dt)")
end   

# Print warnings
function print_warning(s::String)
    printstyled("Warning: ",color=:yellow)
    println(s)    
end

# Print Errors
function print_error(s::String)
    printstyled("Error: ",color=:red)
    println(s)    
end

# Plot results
function PlotResults(;PCHS=nothing,data=nothing,xscale=nothing,graph=nothing,xref=nothing)
    if data===nothing
       print_error("Output dataframe is empty")
    else
        if PCHS===nothing
            print_warning("Empty PCHS")
            PCHS =  build_PHS()
        end
        if xscale===nothing
            xscale = ones(PCHS.nx)
            xl = "x"
        else
            xl = "x * scale" 
            if typeof(xscale)==Matrix{Float64}
                xscale = diag(xscale)
            end
            if length(xscale) != PCHS.nx
                print_warning("Scale is not consistent")
            end
        end
        if graph === nothing
            print_warning("Not selected a graphic. States graphics (X) will plot")
            graph = ["x"]
        end 

        if xref === nothing
            ref = 0
        else
            if length(xref) != PCHS.nx
                print_warning("xref is NOT consistent with number of states")
                ref = 0
            else
                ref = 1
            end
        end 

        # Initialization
        lay = 0        # layout
        plt_var = []   # Data to plot

        for g in graph 
            # If selected H
            if g in ["H","h"]
                plt_H = plot(data.t,data.H,ylabel="H(x)",xlabel="Time (s)",label="")
                lay += 1
                push!(plt_var, plt_H)
            end 
            # If selected X
            if g in ["X","x"]
                plt_X = plot(data.t,data[!,PCHS.xlabel[1]]*xscale[1],ylabel=xl,xlabel="Time (s)",label=PCHS.xlabel[1])
                if PCHS.nx > 1
                    for k = 2:PCHS.nx
                        plt_X = plot!(data.t,data[!,PCHS.xlabel[k]]*xscale[k],label=PCHS.xlabel[k])
                    end
                end

                # Si tiene referencias
                if ref == 1
                    for k = 1:PCHS.nx
                        plt_X = plot!(data.t, fill.(xref[k], length(data.t)), label="x_ref[$k]", linestyle=:dash, color=:gray)
                    end
                end

                lay += 1
                push!(plt_var, plt_X)
            end    
            # If selected U
            if g in ["U","u"]
                plt_U = plot(data.t,data[!,PCHS.ulabel[1]],ylabel="u",xlabel="Time (s)",label=PCHS.ulabel[1])    
                if PCHS.nu > 1
                    for k = 2:PCHS.nu
                        plt_U = plot!(data.t,data[!,PCHS.ulabel[k]],label=PCHS.ulabel[k])
                    end        
                end
                lay += 1
                push!(plt_var, plt_U)
            end  
        end 
    end
    if lay == 0
        print_error("Incorrect data entry")
    else
        plt = plot(plt_var..., layout=(lay, 1))  
        display(plt)
    end 
    return nothing
end

