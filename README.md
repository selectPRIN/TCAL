# TCAL
R codes for Article "Cause of death decomposition of cohort survival comparisons"

These files can be used to calculate TCAL, TCAL by causes of death. 

"ungroup_github.R" is needed as overall age-specific mortality rates (from Human Mortality Database) are defined over single ages groups while deaths by cause are defined over five years age groups.
"ungroup" is a package provided by Pascariu, Rizzi, Schoeley & Danko (see https://cran.r-project.org/web/packages/ungroup/index.html) which allows to disaggregate coarsened deaths by causes so rates and cause-specific number of deaths are defined over the same ages.

"TCAL_causes_functions.R" contains functions that can be used to calculate TCAL (overall or decomposed by cause) on already ungrouped data.

"TCAL_causes_exec.R" is an example of applications of functions in "TCAL_causes_functions.R"
