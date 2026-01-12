using JuMP
using SCIP

# Test Case 1: Simple example that should converge in 2-3 iterations
model = Model(SCIP.Optimizer)

@variable(model, 0 <= x[1:2] <= 1, Bin)

@constraint(model, c1, x[1] + x[2] >= 0.4)
@constraint(model, c2, x[1] - x[2] <= 0.3)

@objective(model, Min, x[1])

write_to_file(model, "testcase/test1.mps")

# Test Case 2: Simpler problem where LP solution is fractional but rounds nicely
model2 = Model(SCIP.Optimizer)

@variable(model2, 0 <= y[1:2] <= 1, Bin)

@constraint(model2, y[1] + y[2] >= 0.9)
@constraint(model2, y[1] - y[2] <= 0.1)

@objective(model2, Min, y[1] + y[2])

write_to_file(model2, "testcase/test2.mps")

# Test Case 3: Taken from feasibility pump literature
model4 = Model(SCIP.Optimizer)

@variable(model, 0 <= z[1:10] <= 1, Bin)

@constraint(model, z[1] + z[2] <= 1.0)
@constraint(model, z[3] + z[4] <= 1.0)
@constraint(model, - z[2] - z[4] + z[5] <= 0.0)
@constraint(model, - z[6] + z[7] <= 0.0)
@constraint(model, - z[6] + z[8] <= 0.0)
@constraint(model, - z[5] - z[7] - z[8] - z[9] <= -1.0)
@constraint(model, - z[8] + z[10] <= 0.0)
@constraint(model, z[9] - z[10] <= 0.0)

@objective(model, Min, - z[1] - z[3] - z[6])

write_to_file(model, "testcase/test3.mps")