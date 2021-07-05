using PkgBenchmark
using NonadiabaticMolecularDynamics

results = benchmarkpkg(NonadiabaticMolecularDynamics)
export_markdown("benchmark.md", results)
