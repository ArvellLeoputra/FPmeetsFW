function addColumn!(display::PumpDisplay, name::String, width::Int, decimals::Int = 0)
    push!(display.column, PumpDisplayColumn(name, width, decimals))
end

function printHeader!(display::PumpDisplay)
    for col in display.column
        print(rpad(col.name, col.width))
    end
    println()
end

function formatValue(val::Float64, decimals::Int)
    if isnan(val)
        return "NaN"
    end

    fixedStr = Printf.format(Printf.Format("%.$(decimals)f"), val)
    dotIdx = findfirst('.', fixedStr)
    intDigits = count(isdigit, fixedStr[1:dotIdx-1])
    if intDigits <= DEF_MAX_INT_DIGITS
        return fixedStr
    end

    # Integer part exceeds the fixed digit cap: fall back to scientific notation,
    # which has a fixed length regardless of magnitude.
    return Printf.format(Printf.Format("%.$(decimals)e"), val)
end

function printRow!(display::PumpDisplay, values...)
    for (col, val) in zip(display.column, values)
        if val isa Float64
            val = formatValue(val, col.decimals)
        end
        print(rpad(string(val), col.width))
    end
    println()
end

# Verbose-mode counterpart to printRow! for a single FW-projection iteration
function printVerboseIteration(
    origObj::Float64,
    projObj::Float64,
    step::Float64,
    nFrac::Int,
    fwIters::Int,
    iterTime::Float64,
    perturbed::Bool,
    restarted::Bool,
    outcome::String
)
    @printf("FW: origObj=%.4f projObj=%.4f step=%.4f nFrac=%d fwIters=%d iterTime=%.4fs P=%s R=%s -> %s\n",
        origObj, projObj, step, nFrac, fwIters, iterTime,
        perturbed ? "*" : " ", restarted ? "*" : " ", outcome)
end