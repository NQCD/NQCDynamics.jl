using FiniteDiff
using NQCModels

function finite_difference_gradient(model::NQCModels.ClassicalModels.ClassicalModel, R)
    f(x) = potential(model, x)
    FiniteDiff.finite_difference_gradient(f, R)
end

function finite_difference_gradient(model::NQCModels.QuantumModels.QuantumModel, R)
    f(x, j, i) = potential(model, x)[j,i]
    grad = [Hermitian(zeros(nstates(model), nstates(model))) for _ in CartesianIndices(R)]
    for i=1:nstates(model)
        for j=1:nstates(model)
            gradient = FiniteDiff.finite_difference_gradient(x->f(x,j,i), R)
            for k in eachindex(R)
                grad[k].data[j,i] = gradient[k]
            end
        end
    end
    grad
end

function test_model(model::NQCModels.Model, atoms; rtol=1e-5) # need a way to ensure that finite_diff and D are always of the same dimension, when by construction this is not always true. Especially if there is only 1 atom
    R = rand(ndofs(model), atoms)
    D = derivative(model, R)
    finite_diff = finite_difference_gradient(model, R)
    if atoms == 1
        return isapprox(finite_diff[1], D, rtol=rtol)
    end
    return isapprox(finite_diff, D, rtol=rtol)
end

function test_model(model::NQCModels.QuantumModels.GatesHollowayElbow, atoms; rtol=1e-5)
    R = rand(ndofs(model), atoms)
    D = derivative(model, R)
    finite_diff = finite_difference_gradient(model, R)
    return isapprox(finite_diff, D, rtol=rtol)
end

function test_model(model::NQCModels.WideBandBath, atoms; rtol=1e-5) # need a way to ensure that finite_diff and D are always of the same dimension, when by construction this is not always true. Especially if there is only 1 atom
    R = rand(ndofs(model), atoms)
    D = derivative(model, R)
    finite_diff = finite_difference_gradient(model, R)
    return isapprox(finite_diff, D, rtol=rtol)
end

function test_model(model::NQCModels.QuantumModels.AndersonHolstein, atoms; rtol=1e-5)
    R = rand(ndofs(model), atoms)
    D = derivative(model, R)
    finite_diff = finite_difference_gradient(model, R)
    return isapprox(finite_diff, D, rtol=rtol)
end

function test_model(model::NQCModels.QuantumModels.DoubleWell, atoms; rtol=1e-5)
    R = rand(ndofs(model), atoms)
    D = derivative(model, R)
    finite_diff = finite_difference_gradient(model, R)
    return isapprox(finite_diff, D, rtol=rtol)
end

function test_model(model::NQCModels.QuantumModels.Scattering1D, atoms; rtol=1e-5)
    R = rand(ndofs(model), atoms)
    D = derivative(model, R)
    finite_diff = finite_difference_gradient(model, R)
    return isapprox(finite_diff, D, rtol=rtol)
end

function test_model(model::NQCModels.QuantumModels.AdiabaticStateSelector, atoms; rtol=1e-5)
    R = rand(ndofs(model), atoms)
    D = derivative(model, R)
    finite_diff = finite_difference_gradient(model, R)
    return isapprox.(finite_diff, D, rtol=rtol) |> any
end

function test_model(model::NQCModels.FrictionModels.ClassicalFrictionModel, atoms)
    R = rand(ndofs(model), atoms)
    D = derivative(model, R)
    friction(model, R)
    return finite_difference_gradient(model, R) ≈ D
end

function test_model(model::NQCModels.ClassicalModels.CSVModel_1D, atoms; rtol=1e-5)
    R = rand(ndofs(model), atoms)
    D = derivative(model, R)
    finite_diff = finite_difference_gradient(model, R)
    return isapprox(finite_diff, D, rtol=rtol)
end