using DataFrames, LaTeXStrings

function bodeplot_SL(G)
    zs_data, ps_data, k_data  = zpkdata(G)
    zs, ps, k = first(zs_data), first(ps_data), first(k_data)

    PZ  = DataFrame(value = Float64[], type = Symbol[])
    for p in ps
        if (imag(p) ≉ 0.0)
            error("Found a compex pole: $p. Does not handle complex poles or zeros")
        end
        # Assume that the system is stable
        p = abs(p)
        push!(PZ, (value = log10(p) - 1, type = :p_minus))
        push!(PZ, (value = log10(p)    , type = :p))
        push!(PZ, (value = log10(p) + 1, type = :p_plus))
    end 
    for z in zs
        if (imag(z) ≉ 0.0)
            error("Found a compex zero: $z. Does not handle complex poles or zeros")
        end
        z = abs(z)
        push!(PZ, (value = log10(z) - 1, type = :z_minus))
        push!(PZ, (value = log10(z)    , type = :z))
        push!(PZ, (value = log10(z) + 1, type = :z_plus))
    end 


    sort!(PZ, :value)

    df = DataFrame(log10ω = Float64[], gain = Float64[], phase = Float64[])

    gain_slope = 0
    phase_slope = 0
    phase = 0

    while(PZ[1, :value] == -Inf)
        if (PZ[1, :type] == :p)
            gain_slope -= 20
            phase -= 90
        elseif(PZ[1, :type] == :z)
            gain_slope += 20
            phase += 90
        end
        popfirst!(PZ)
    end

    # Add points at the beginning and end
    pushfirst!(PZ, (value = PZ[1, :value] - 1, type = :first))
    push!(PZ, (value = PZ[end, :value] + 1, type = :last))


    k *= prod(abs, filter(x -> x ≉ 0.0, zs))
    k /= prod(abs, filter(x -> x ≉ 0.0, ps))

    log10ω = PZ[1, :value]
    gain = 20*log10(k) + log10ω * gain_slope 

    push!(df, (log10ω=log10ω, gain=gain, phase=phase))

    for (pz, type) in eachrow(PZ)
        current_log10ω = log10ω
        current_gain   = gain
        current_phase  = phase

        log10ω = pz
        gain   = current_gain  + gain_slope  * (log10ω - current_log10ω)
        phase  = current_phase + phase_slope * (log10ω - current_log10ω)
        
        push!(df, (log10ω=log10ω, gain=gain, phase=phase))
        
        if (type == :z)
            gain_slope += 20
        elseif (type == :p)
            gain_slope -= 20
        elseif (type == :z_minus || type == :p_plus)
            phase_slope += 45
        elseif (type == :z_plus || type == :p_minus)
            phase_slope -= 45
        end 
    end

    return combine(groupby(df, :log10ω),
                   :gain => last, :phase => last,
                   renamecols = false)
end

function get_yticks(df, col, step, suffix)
    data_max = maximum(df[!, col])
    data_min = minimum(df[!, col])
    yticks = range(floor(Int, data_min/step)*step,
                   ceil(Int, data_max/step)*step,
                   step=step) |> collect
    ylabels = string.(yticks, suffix)

    return (yticks, ylabels)
end
    
function get_xticks(df)
    col = :log10ω
    step = 1
    data_max = maximum(df[!, col])
    data_min = minimum(df[!, col])
    xticks = range(floor(Int, data_min/step)*step,
                   ceil(Int, data_max/step)*step,
                   step=step) |> collect
    xlabels = [L"10^{%$xtick}" for xtick in xticks]

    return (xticks, xlabels)
end
    
