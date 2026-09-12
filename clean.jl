# Discards a JIT-warm-up solve so the second solve's timing is trustworthy.
using FPmeetsFW

fileName, config, resultsDir = loadConfig(ARGS)
instanceName = basename(fileName)

# Warm-up: same instance, same config, discarded - no results written
println("Warming up on $instanceName...")
redirect_stdout(devnull) do
    runInstance(fileName, config, "")
end

# Real run: fully compiled now, this timing is trustworthy
println("Running $instanceName...")
runInstance(fileName, config, resultsDir)
