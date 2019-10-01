module vapor_pressure
cd(dirname(@__FILE__))
include("./chemical.jl")
include("./thermo_model.jl")
using .chemical: Chemical
using .thermo_model: TDependentModelHandler, TDependentModel
using DataFrames: DataFrame, names
using CSV: read

function sample_dict(name, symbol; exclude)
    table = read(name)
    oldcols = names(table)
    mask = (oldcols.!=symbol)
    other_masks = [oldcols .!= i for i in exclude]
    mask = 
    oldcols = oldcols[ .== (oldcols.!= :Name) .== (oldcols.!= :Chemical)]
    data = hcat((table[!, col] for col in oldcols)...)
    newcols = [replace(i, r"-" => "") for i in table[!, symbol]]
    Dict(col => data[i, :] for (i, col) in enumerate(newcols))
end

cd("./Vapor Pressure")
const name = "vapor pressure"
const WagnerMcGarry = sample_dict(read("Wagner Original McGarry.tsv"), :CASRN, exclude=[:Name])
const Wagner = sample_dict(read("Wagner Collection Poling.tsv"), :CASRN, exclude=[:Name])
const Antoine = sample_dict(read("Antoine Collection Poling.tsv"), :CASRN, exclude=[:Chemical])
const AntoineExtended = sample_dict(read("Antoine Extended Collection Poling.tsv"), :CASRN, exclude=[:Chemical])
const Perrys2_8 = sample_dict(read("Table 2-8 Vapor Pressure of Inorganic and Organic Liquids.tsv"), :CAS, exclude=[:Chemical])
const VDI_PPDS_3 = sample_dict(read("VDI PPDS Boiling temperatures at different pressures.tsv"), :CAS, exclude=[:Chemical])
const WagnerMcGarry_keys = keys(WagnerMcGarry)
const Antoine_keys = keys(Antoine)
const Wagner_keys = keys(Wagner)
const AntoineExtended_keys = keys(AntoineExtended)
const Perrys2_8_keys = keys(Perrys2_8)
const VDI_PPDS_3 = keys(VDI_PPDS_3_keys)

function AntoineModel(A::Float64, B::Float64, C::Float64, Tmin::Float64, Tmax::Float64)
    TDependentModel("Antoine", T::Float64 -> 10.0^(A - B / (T + C)), Tmin, Tmax)
end

function AntoineExtendedModel(Tc::Float64, to::Float64, A::Float64, B::Float64, C::Float64, n::Float64, E::Float64, F::Float64, Tmin::Float64, Tmax::Float64)
    function vapor_pressure(T::Float64)
        x = max((T - to - 273.15) / Tc, 0.0)
        10.0^(A - B / (T + C) + 0.43429 * x^n + E * x^8 + F * x^12)
    end
    TDependentModel("Antoine extended", vapor_pressure, Tmin, Tmax)
end

function WagnerMcGarryModel(A::Float64, B::Float64, C::Float64, D::Float64, Tc::Float64, Pc::Float64, Tmin::Float64, Tmax::Float64)
    function vapor_pressure(T::Float64)
        Tr = T / Tc
        tau = 1.0 - Tr
        Pc * exp((A * tau + B * tau^1.5 + C * tau^3 + D * tau^6) / Tr)
    end
    TDependentModel("Wagner McGarry", vapor_pressure, Tmin, Tmax)
end

function WagnerModel(Tc::Float64, Pc::Float64, A::Float64, B::Float64, C::Float64, D::Float64, Tmin::Float64, Tmax::Float64)
    function vapor_pressure(T::Float64)
        Tr = T / Tc
        τ = 1.0 - Tr
        Pc * exp((A*τ + B*τ^1.5 + C*τ^2.5 + D*τ^5) / Tr)
    end
    TDependentModel("Wagner", vapor_pressure, Tmin, Tmax)
end

function BoilingCriticalRelationModel(Tb::Float64, Tc::Float64, Pc::Float64, Tmin::Float64, Tmax::Float64)
    Tbr = Tb / Tc
    h = Tbr * log(Pc / 101325.0) / (1 - Tbr)
    TDependentModel("Boiling-critical relation", T::Float64 -> exp(h * (1 - Tc / T)) * Pc, Tmin, Tmax)
end

function LeeKeslerModel(Tc::Float64, Pc::Float64, ω::Float64, Tmin::Float64, Tmax::Float64)
    function vapor_pressure(T::Float64)
        Tr = T / Tc
        Tra = Tr^6
        logTr = log(Tr)
        f0 = 5.92714 - 6.09648 / Tr - 1.28862 * logTr + 0.169347 * Tra
        f1 = 15.2518 - 15.6875 / Tr - 13.4721 * logTr + 0.43577 * Tra
        exp(f0 + ω * f1) * Pc
    end
    TDependentModel("Lee Kesler", vapor_pressure, Tmin, Tmax)
end

function AmbroseWaltonModel(Tc::Float64, Pc::Float64, ω::Float64, Tmin::Float64, Tmax::Float64)
    function vapor_pressure(T::Float64)
        Tr = T / Tc
        τ = 1 - Tr
        τa = τ^1.5
        τb = τ^2.5
        τc = τ^5
        f0 = -5.97616 * τ + 1.29874 * τa - 0.60394 * τb - 1.06841 * τc
        f1 = -5.03365 * τ + 1.11505 * τa - 5.41217 * τb - 7.46628 * τc
        f2 = -0.64771 * τ + 2.41539 * τa - 4.26979 * τb + 3.25259 * τc
        Pc * exp((f0 + f1 * ω + f2 * ω^2) / Tr)
    end
    TDependentModel("Ambrose Walton", vapor_pressure, Tmin, Tmax)
end

function SanjariModel(Tc::Float64, Pc::Float64, ω::Float64, Tmin::Float64, Tmax::Float64)
    function vapor_pressure(T::Float64)
        Tr = T / Tc
        logTr = log(Tr)
        Ta = Tr^1.9
        f0 = 6.83377 + -5.76051 / Tr + 0.90654 * logTr + -1.16906 * Ta
        f1 = 5.32034 + -28.1460 / Tr + -58.0352 * logTr + 23.57466 * Ta
        f2 = 18.19967 + 16.33839 / Tr + 65.6995 * logTr + -35.9739 * Ta
        Pc * exp(f0 + f1 * ω + f2 * ω^2)
    end
    TDependentModel("Sanjari", vapor_pressure, Tmin, Tmax)
end

function EOSModel(EOS, Tmin, Tmax)
    TDependentModel("EOS", T -> EOS.Psat(T), Tmin, Tmax)
end

function EdalatModel(Tc::Float64, Pc::Float64, ω::Float64, Tmin::Float64, Tmax::Float64)
    function vapor_pressure(T::Float64)
        τ = 1.0 - T / Tc
        a = -6.1559 - 4.0855 * ω
        c = -0.8747 - 7.8874 * ω
        d = 1.0 / (-0.4893 - 0.9912 * ω + 3.1551 * ω^2)
        b = 1.5737 - 1.0540 * ω - 4.4365E-3 * d
        lnPr = (a * τ + b * τ^1.5 + c * τ^3.0 + d * τ^6.0) / (1.0 - τ)
        exp(lnPr) * Pc
    end
    TDependentModel("Edalat", vapor_pressure, Tmin, Tmax)
end

function load_vapor_pressure_models(chemical::Chemical)
    chemical.vapor_pressure = vapor_pressure_model_handler(chemical.CAS, chemical.Tc, chemical.Pc, chemical.ω)
end

function vapor_pressure_model_handler(CAS; Tb = nothing, Tc = nothing, Pc = nothing, ω = nothing, eos = nothing, models = nothing)
    if models == nothing
        models = Vector{TDependentModel}()
    end
    if CAS in WagnerMcGarry_keys
        A, B, C, D, Pc, Tc, Tmin = WagnerMcGarry[CAS]
        Tmax = Tc
        push!(models, WagnerMcGarryModel(A, B, C, D, Tc, Pc, Tmin, Tmax))
    end
    if CAS in Wagner_keys
        A, B, C, D, Tc, Pc, Tmin, Tmax = WagnerPoling[CAS]
        # Some Tmin values are missing; Arbitrary choice of 0.1 lower limit
        if isnan(Tmin)
            Tmin = Tmax * 0.1
        end
        push!(models, WagnerModel(A, B, C, D, Tc, Pc, Tmin, Tmax))
    end
    if CAS in AntoineExtended_keys
        A, B, C, Tc, to, n, E, F, Tmin, Tmax = AntoineExtended[CAS]
        push!(models, AntoineExtendedModel(A, B, C, Tc, to, n, E, F, Tmin, Tmax))
    end
    if CAS in Antoine_keys
        A, B, C, Tmin, Tmax = Antoine[CAS]
        push!(models, AntoineModel(A, B, C, Tmin, Tmax))
    end
    if CAS in Perrys2_8_keys
        # push!(models, DIPPR_PERRY_8E)
        # C1, C2, C3, C4, C5, Tmin, Tmax = Perrys2_8[!, CAS]
    end
    if CAS in VDI_PPDS_3_keys
        Tm, Tc, Pc, A, B, C, D = VDI_PPDS_3[CAS]
        push!(models, WagnerModel(Tc, Pc, A, B, C, D, Tmin, Tmax))
    end
    if all([Tb, Pc, Pc] .!= nothing)
        push!(models, BoilingCriticalRelationModel(Tb, Tc, Pc))
    end
    if all([Tc, Pc, ω] .!= nothing)
        push!(models, LeeKeslerModel(Tc, Pc, ω))
        push!(models, AmbroseWaltonModel(Tc, Pc, ω))
        push!(models, SanjariModel(Tc, Pc, ω))
        push!(models, EdalatModel(Tc, Pc, ω))
    end
    if eos != nothing
        Tmin = 0.01
        Tmax = Tc
        if eos
            push!(models, EOSModel(eos, Tmin, Tmax))
        end
    end
    if length(models) == 0
        throw(DomainError(CAS, "no vapor pressure models available"))
    end
    acting_model = models[1]
    TDependentModelHandler(name, models, acting_model)
end



end
