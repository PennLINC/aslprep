# AGENTS.md -- ASLPrep

This file provides instructions for AI coding agents and human maintainers working on **ASLPrep**, a BIDS App for preprocessing Arterial Spin Labeling (ASL) MRI data.

---

## Shared Instructions (All PennLINC BIDS Apps)

The following conventions apply equally to **qsiprep**, **qsirecon**, **xcp_d**, and **aslprep**. All four are PennLINC BIDS Apps built on the NiPreps stack.

### Ecosystem Context

- These projects belong to the [NiPreps](https://www.nipreps.org/) ecosystem and follow its community guidelines.
- Core dependencies include **nipype** (workflow engine), **niworkflows** (reusable workflow components), **nireports** (visual reports), **pybids** (BIDS dataset querying), and **nibabel** (neuroimaging I/O).
- All four apps are containerized via Docker and distributed on Docker Hub under the `pennlinc/` namespace.
- Contributions follow the [NiPreps contributing guidelines](https://www.nipreps.org/community/CONTRIBUTING/).

### Architecture Overview

Every PennLINC BIDS App follows this execution flow:

```
CLI (parser.py / run.py)
  -> config singleton (config.py, serialized as ToML)
    -> workflow graph construction (workflows/*.py)
      -> Nipype interfaces (interfaces/*.py)
        -> BIDS-compliant derivative outputs
```

- **CLI** (`<pkg>/cli/`): `parser.py` defines argparse arguments; `run.py` is the entry point; `workflow.py` builds the execution graph; `version.py` handles `--version`.
- **Config** (`<pkg>/config.py`): A singleton module with class-based sections (`environment`, `execution`, `workflow`, `nipype`, `seeds`). Config is serialized to ToML and passed across processes via the filesystem. Access settings as `config.section.setting`.
- **Workflows** (`<pkg>/workflows/`): Built using `nipype.pipeline.engine` (`pe.Workflow`, `pe.Node`, `pe.MapNode`). Use `LiterateWorkflow` from `niworkflows.engine.workflows` for auto-documentation. Every workflow factory function must be named `init_<descriptive_name>_wf`.
- **Interfaces** (`<pkg>/interfaces/`): Custom Nipype interfaces wrapping external tools or Python logic. Follow standard Nipype patterns: define `_InputSpec` / `_OutputSpec` with `BaseInterfaceInputSpec` / `TraitedSpec`, implement `_run_interface()`.
- **Utilities** (`<pkg>/utils/`): Shared helper functions. BIDS-specific helpers live in `utils/bids.py`.
- **Reports** (`<pkg>/reports/`): HTML report generation using nireports.
- **Data** (`<pkg>/data/`): Static package data (config files, templates, atlases). Accessed via `importlib.resources` or the `acres` package.
- **Tests** (`<pkg>/tests/`): Pytest-based. Unit tests run without external data. Integration tests are gated behind pytest markers and are skipped by default.

### Workflow Authoring Rules

1. Every workflow factory function must be named `init_<name>_wf` and return a `Workflow` object.
2. Use `LiterateWorkflow` (from `niworkflows.engine.workflows`) to enable automatic workflow graph documentation.
3. Define `inputnode` and `outputnode` as `niu.IdentityInterface` nodes to declare the workflow's external API.
4. Connect nodes using `workflow.connect([(source, dest, [('out_field', 'in_field')])])` syntax.
5. Add `# fmt:skip` after multi-line `workflow.connect()` calls to prevent ruff from reformatting them.
6. Include a docstring with `Workflow Graph` and `.. workflow::` Sphinx directive for auto-generated documentation.
7. Use `config` module values (not function parameters) for global settings inside workflow builders.

### Interface Conventions

1. Input/output specs use Nipype traits (`File`, `traits.Bool`, `traits.Int`, etc.).
2. `mandatory = True` for required inputs; provide `desc=` for all traits.
3. Implement `_run_interface(self, runtime)` -- never `run()`.
4. Return `runtime` from `_run_interface`.
5. Set outputs via `self._results['field'] = value`.

### Config Module Usage

```python
from <pkg> import config

# Read a setting
work_dir = config.execution.work_dir

# Serialize to disk
config.to_filename(path)

# Load from disk (in a subprocess)
config.load(path)
```

The config module is the single source of truth for runtime parameters. Never pass global settings as function arguments when they are available via config.

### Testing Conventions

- **Unit tests**: Files named `test_*.py` in `<pkg>/tests/`. Must not require external neuroimaging data or network access.
- **Integration tests**: Decorated with `@pytest.mark.<marker_name>`. Excluded by default via `addopts` in `pyproject.toml`. Require Docker or pre-downloaded test datasets.
- **Fixtures**: Defined in `conftest.py`. Common fixtures include `data_dir`, `working_dir`, `output_dir`, and `datasets`.
- **Coverage**: Configured in `pyproject.toml` under `[tool.coverage.run]` and `[tool.coverage.report]`.

### Documentation

- Built with Sphinx using `sphinx_rtd_theme`.
- Source files in `docs/`.
- Workflow graphs are auto-rendered via `.. workflow::` directives that call `init_*_wf` functions.
- API docs generated via `sphinxcontrib-apidoc`.
- Bibliography managed with `sphinxcontrib-bibtex` and `boilerplate.bib`.

### Docker

- Each app has a base image with runtime dependencies and a main Dockerfile that installs the Python environment and the app itself.
- ASLPrep and fMRIPrep have migrated to **Pixi**-based multi-stage Docker builds (base has no Python; Pixi installs conda + PyPI deps from `pixi.lock`). The other repos still use micromamba + `pip install`.
- Base image naming: ASLPrep uses `pennlinc/aslprep-base:<YYYYMMDD>` (date-based); other repos use `pennlinc/<pkg>_build:<version>`.
- Labels follow the `org.label-schema` convention.

### Release Process

- Versions are derived from git tags via `hatch-vcs` (VCS-based versioning).
- GitHub Releases use auto-generated changelogs configured in `.github/release.yml`.
- Release categories: Breaking Changes, New Features, Deprecations, Bug Fixes, Other.
- Docker images are built and pushed via CI on tagged releases.

### Code Style

- **Formatter**: `ruff format` (target: all four repos).
- **Linter**: `ruff check` with an extended rule set (F, E, W, I, UP, YTT, S, BLE, B, A, C4, DTZ, T10, EXE, FA, ISC, ICN, PT, Q).
- **Import sorting**: Handled by ruff's `I` rule (isort-compatible).
- **Pre-commit**: Uses `ruff-pre-commit` hooks for both linting and formatting.
- **Black is disabled**: `[tool.black] exclude = ".*"` in repos that have migrated to ruff.

### BIDS Compliance

- All outputs must conform to the [BIDS Derivatives](https://bids-specification.readthedocs.io/en/stable/derivatives/introduction.html) specification.
- Use `pybids.BIDSLayout` for querying input datasets.
- Use `DerivativesDataSink` (from the project's interfaces or niworkflows) for writing BIDS-compliant output files.
- Entity names, suffixes, and extensions must match the BIDS specification.

---

## ASLPrep-Specific Instructions

### Project Overview

ASLPrep is a BIDS App for preprocessing Arterial Spin Labeling (ASL) perfusion MRI data. It handles:
- Anatomical preprocessing (via smriprep/fmriprep infrastructure)
- ASL reference image estimation
- Head motion correction for ASL time series
- Susceptibility distortion correction (via sdcflows)
- Cerebral blood flow (CBF) computation for single-PLD and multi-PLD acquisitions
- Support for PCASL and PASL labeling schemes
- Confound estimation
- Parcellation-based CBF summarization
- QC metric aggregation (`aslprep-combineqc` CLI tool)

### Repository Details

| Item | Value |
|------|-------|
| Package name | `aslprep` |
| Default branch | `main` |
| Entry point | `aslprep.cli.run:main` |
| Python requirement | `>=3.10` |
| Build backend | hatchling + hatch-vcs |
| Linter | ruff ~= 0.15.0 |
| Pre-commit | Yes (ruff v0.6.2) |
| Tox | Yes |
| Docker base | `pennlinc/aslprep-base:<YYYYMMDD>` |
| Dockerfile | Multi-stage pixi build |

### Key Directories

- `aslprep/workflows/asl/`: ASL-specific workflows organized by processing stage:
  - `fit.py`: Model fitting
  - `base.py`: Main ASL workflow orchestration
  - `cbf.py`: Cerebral blood flow computation
  - `confounds.py`: Confound estimation
  - `hmc.py`: Head motion correction
  - `reference.py`: Reference image estimation
  - `apply.py`: Transform application
  - `resampling.py`: Resampling to standard spaces
  - `outputs.py`: BIDS derivative output writing
  - `plotting.py`: QC visualization
- `aslprep/workflows/base.py`: Top-level workflow (subject iteration)
- `aslprep/workflows/segmentation.py`: Brain segmentation
- `aslprep/interfaces/cbf.py`: CBF computation interfaces
- `aslprep/interfaces/parcellation.py`: Atlas parcellation interfaces
- `aslprep/utils/cbf.py`: CBF calculation utilities
- `aslprep/utils/asl.py`: ASL-specific helper functions

### Heavy Dependency on fMRIPrep/NiPreps Stack

ASLPrep has the deepest dependency chain of the four repos:
- `fmriprep ~= 25.2.2`: Reuses fMRIPrep's anatomical and fieldmap workflows
- `sdcflows >= 2.15.0`: Susceptibility distortion correction
- `smriprep >= 0.19.2`: Structural MRI preprocessing
- `niworkflows >= 1.14.4`: Shared workflow components
- `nitransforms >= 25.0.1`: Spatial transforms

When updating these dependencies, check for breaking API changes across the stack.

### Version Management

ASLPrep imports `__version__` directly from `_version.py`.

### Test Extras Naming

ASLPrep names its test dependencies `tests` (plural) in `pyproject.toml`, consistent with the other three repos.

### Template Fetching

ASLPrep's tox config references `python scripts/fetch_templates.py` as a `commands_pre` step before tests. This script downloads the required TemplateFlow templates and is also used by the Docker build's `templates` stage.

### Linting Notes

ASLPrep has 4 suppressed ruff rules:
- `UP031`: Old-style string formatting (to be fixed incrementally)
- `S311`: Random not for crypto (intentional)
- `ISC001`: Implicit string concatenation (conflicts with formatter)
- `S603`: Subprocess with shell=True (trusted commands only)

### ASL-Specific Terminology

When working on ASL workflows, know these domain terms:
- **PLD** (Post-Labeling Delay): Time between labeling and image acquisition
- **CBF** (Cerebral Blood Flow): The primary output metric
- **PCASL** (Pseudo-Continuous ASL): A labeling technique
- **PASL** (Pulsed ASL): Another labeling technique
- **Single-PLD vs Multi-PLD**: Whether one or multiple delay times are used
- **Label-control pairs**: ASL acquires alternating labeled and control images

---

## Cross-Project Development Roadmap

This roadmap covers harmonization work across all four PennLINC BIDS Apps (qsiprep, qsirecon, xcp_d, aslprep) to reduce maintenance burden.

### Phase 1: Bring qsirecon to parity

1. **Migrate qsirecon from flake8+black+isort to ruff** -- copy the `[tool.ruff]` config from xcp_d's `pyproject.toml` and remove `[tool.black]`, `[tool.isort]`, `[tool.flake8]` sections.
2. **Add `.pre-commit-config.yaml` to qsirecon** -- identical to the config used by qsiprep, xcp_d, and aslprep.
3. **Add `tox.ini` to qsirecon** -- copy from qsiprep or xcp_d (they are identical).
4. **Add `.github/dependabot.yml` to qsirecon**.
5. **Reformat qsirecon codebase** -- run `ruff format` to switch from double quotes to single quotes.

### Phase 2: Standardize across all four repos

6. **Rename qsiprep default branch** from `master` to `main` and update `.github/workflows/lint.yml`.
7. ~~**Rename aslprep test extras** from `test` to `tests`~~ -- **Done**.
8. **Converge on version management** -- recommend the simpler `_version.py` direct-import pattern (used by qsiprep/qsirecon). Migrate xcp_d and aslprep away from `__about__.py`.
9. **Pin the same ruff version** in all four repos' dev dependencies and `.pre-commit-config.yaml`.
10. **Harmonize ruff ignore lists** -- adopt xcp_d's minimal set (`S105`, `S311`, `S603`) as the target; fix suppressed rules in qsiprep and aslprep incrementally.

### Phase 3: Shared infrastructure

11. **Extract a reusable GitHub Actions workflow** for lint + codespell + build checks, hosted in a shared repo (e.g., `PennLINC/.github`).
12. **Standardize Dockerfile patterns** -- ASLPrep and fMRIPrep have migrated to Pixi-based multi-stage builds. Migrate qsiprep, qsirecon, and xcp_d to the same pattern.
13. **Create a shared `pennlinc-style` package or cookiecutter template** providing `pyproject.toml` lint/test config, `.pre-commit-config.yaml`, `tox.ini`, and CI workflows.
14. **Evaluate `nipreps-versions` calver** -- the `raw-options = { version_scheme = "nipreps-calver" }` line is commented out in all four repos. Decide whether to adopt it.

