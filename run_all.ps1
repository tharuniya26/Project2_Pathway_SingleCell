$ErrorActionPreference = "Stop"

$ProjectRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$Python = Join-Path $ProjectRoot ".venv\Scripts\python.exe"

if (-not (Test-Path -LiteralPath $Python)) {
    throw "Python environment not found at $Python. See README.md for setup."
}

Push-Location $ProjectRoot
try {
    & $Python "src\run_adjusted_differential_expression.py" --n-cpus 1 --quiet
    & $Python "src\run_comparative_analysis.py"
    & $Python "src\validate_project.py"
}
finally {
    Pop-Location
}
