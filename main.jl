# Load FPmeetsFW as the active-project package (run with --project=<this dir>) so
# a prebuilt sysimage's compiled code is actually reused. See scripts/sysimage/.
using FPmeetsFW

fileName, config, resultsDir = loadConfig(ARGS)
runInstance(fileName, config, resultsDir)
