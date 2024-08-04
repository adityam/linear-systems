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

function gain_crossover(df)
    initial_gain = df[1, :gain]
    compare = (initial_gain ≥ 0) ? (<) : (≥)

    done = false
    ω_lower = df[1, :log10ω]
    g_lower = df[1, :gain]

    ω_upper = df[end, :log10ω]
    g_upper = df[end, :gain]
    for (log10ω, gain, phase) = eachrow(df)
        if compare(gain,0)
            ω_upper = log10ω
            g_upper = gain
            done = true
            break
        end
        ω_lower = log10ω
        g_lower = gain
    end

    if (!done)
        @debug("0dB line does not intersect in the range of dataframe")
        return nothing
    end

    slope = (g_lower - g_upper)/(ω_lower - ω_upper)
    ω = ω_upper - g_upper/slope
end

function phase_crossover(df)
    initial_phase = df[1, :phase]
    compare = (initial_phase ≥ -180) ? (≤) : (≥)

    done = false
    ω_lower = df[1, :log10ω]
    ϕ_lower = df[1, :phase]

    ω_upper = df[end, :log10ω]
    ϕ_upper = df[end, :phase]
    for (log10ω, gain, phase) = eachrow(df)
        if compare(phase,-180)
            ω_upper = log10ω
            ϕ_upper = phase
            done = true
            break
        end
        ω_lower = log10ω
        ϕ_lower = phase
    end

    if (!done)
        @debug("-180° line does not intersect in the range of dataframe")
        return nothing
    end

    slope = (ϕ_lower - ϕ_upper)/(ω_lower - ω_upper)
    ω = ω_upper + (-180 - ϕ_upper)/slope
end

function get_value(df, ω)
    result = searchsorted(df.log10ω, ω)
    
    idx_lower = last(result)
    idx_upper = first(result)

    if (idx_lower == 0 || idx_upper > nrow(df))
        @debug("ω not in range")
        return (nothing, nothing)
    end

    (ω_lower, g_lower, ϕ_lower) = df[idx_lower, :]
    (ω_upper, g_upper, ϕ_upper) = df[idx_upper, :]

    if (idx_lower == idx_upper)
        return (g_lower, ϕ_lower)
    end

    g_slope = (g_lower - g_upper)/(ω_lower - ω_upper)
    ϕ_slope = (ϕ_lower - ϕ_upper)/(ω_lower - ω_upper)

    g = g_lower - g_slope*(ω_lower - ω) 
    ϕ = ϕ_lower - ϕ_slope*(ω_lower - ω) 
    return (g, ϕ)
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
    
