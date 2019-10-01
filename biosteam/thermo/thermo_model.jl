module thermo_model
abstract type ThermoModel end
abstract type ThermoModelHandler end
import Base.show
import Base.typeof

function typename(obj)
    type = typeof(obj)
    split(string(type), ".")[end]
end

show(io::IO, model::Union{ThermoModelHandler, ThermoModel}) = print(io, "<$(typename(model)): $(uppercasefirst(model.name))>")

struct TDependentModel <: ThermoModel
    name::String
    calculate::Function
    Tmin::Float64
    Tmax::Float64
end

struct TPDependentModel <: ThermoModel
    name::String
    calculate::Function
    Tmin::Float64
    Tmax::Float64
    Pmin::Float64
    Pmax::Float64
end

struct TDependentModelHandler <: ThermoModelHandler
    name::String
    models::Array{TDependentModel, 1}
    acting_model::ThermoModel
end

struct TPDependentModelHandler <: ThermoModelHandler
    name::String
    models::Array{TPDependentModel, 1}
    acting_model::ThermoModel
end

is_within_domain(model::TDependentModel, T::Float64) = model.Tmin < T < model.Tmax
is_within_domain(model::TPDependentModel, T::Float64, P::Float64) = (model.Tmin < T < model.Tmax) & (model.Pmin < P < model.Pmax)

function (handler::TDependentModelHandler)(T::Float64)
    model = handler.acting_model
    if is_within_domain(model, T)
        return model.calculate(T)
    else
        other_models = handler.models[[model==i for i in handler.models]]
        for model in other_models
            if is_within_domain(model, T)
                handler.acting_model = model
                return model.calculate(T, P)
            end
        end
        throw(DomainError(T, "no $(handler.name) model at requested T"))
    end
end

function (handler::TPDependentModelHandler)(T::Float64, P::Float64)
    model = handler.acting_model
    if is_within_domain(model, T, P)
        return model.calculate(T, P)
    else
        other_models = handler.models[handler.models .!= model]
        for model in other_models
            if is_within_domain(model, T, P)
                handler.acting_model = model
                return model.calculate(T, P)
            end
        end
        throw(DomainError((T, P), "no $(handler.name) model at requested T and P"))

    end
end



end
