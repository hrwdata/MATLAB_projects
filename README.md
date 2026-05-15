# MATLAB Projects

## Purpose

This repository collects several small MATLAB research and teaching projects covering ordinary differential equations, numerical dynamics, optimization, and simple epidemiological modeling.

## Repository Layout

- `src/`: reusable MATLAB functions grouped by topic
- `examples/`: runnable demos and Live Scripts
- `tests/`: `matlab.unittest` regression and smoke tests
- `docs/`: public documentation for setup, algorithms, testing, and data provenance
- `sample-data/`: small example datasets suitable for redistribution
- `archive/`: temporary preservation area for migration artifacts and generated exports under review
- `projects/`: preserved original snapshot layout kept during the first reorganization pass

## Algorithms Included

- Euler's method for initial value problems
- Direction-field and bifurcation plotting helpers
- Armijo, Goldstein, and Wolfe line search methods
- Partial conjugate gradient optimization
- Bounded simplex and primal-dual linear programming routines
- SIR model fitting with weekly aggregation of case data

## Requirements

- MATLAB with ODE solvers for most examples
- Symbolic Math Toolbox for `pcgrad` and `sir_fit`
- A valid MATLAB license is required to run automated tests in MATLAB

## Installation / Setup

1. Clone the repository.
2. In MATLAB, add the reusable source tree to the path:

```matlab
addpath(genpath('src'));
addpath(genpath('examples'));
addpath(genpath('sample-data'));
```

3. Review [docs/setup.md](docs/setup.md) for environment notes and project layout guidance.

## Minimal Reproducible Example

```matlab
repoRoot = pwd;
addpath(genpath(fullfile(repoRoot, 'src')));
originalDir = pwd;
cleanupObj = onCleanup(@() cd(originalDir));
cd(fullfile(repoRoot, 'src', 'numerics'));
[x, y, f, results_table] = eulers(@(x,y) x + y, 0, 1, 1, 0.1);
disp(results_table(1:3,:))
```

## Running Examples

- Open the Live Scripts under `examples/tutorial/`, `examples/bifurcation/`, and `examples/epidemiology/` in MATLAB.
- Run the MATLAB example scripts in `examples/optimization/` from the repository root after adding `src/` to the path.

## Running Tests

```matlab
repoRoot = pwd;
addpath(genpath(fullfile(repoRoot, 'src')));
originalDir = pwd;
cleanupObj = onCleanup(@() cd(originalDir));
cd(fullfile(repoRoot, 'tests'));
results = runtests('.');
table(results)
```

See [docs/testing.md](docs/testing.md) for test scope and toolbox notes.

## Data Provenance

The repository includes only a small sample COVID dataset intended for examples. See [docs/data-provenance.md](docs/data-provenance.md) for details.

## Citation

If you use this repository in teaching or research, cite the repository URL and the specific commit or release used.

## License

No license file has been added in this pass. Confirm the intended open-source license before publishing or accepting external contributions.
