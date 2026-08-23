---
title: Current System Architecture
status: CURRENT_VERIFIED
owner: Andrea Zedda
audience: engineers, researchers and reviewers
last_verified_at: 2026-08-24
last_verified_git_sha: bf097810b337dc6b766cda04497005670cd96513
source_of_truth: repository tree, Django configuration, model registry and current contracts
---

# Current system architecture

> **CURRENT architecture.** This page describes the repository at the verified SHA. It is not the target virtual-laboratory architecture.

## System context

```mermaid
flowchart LR
    U["Research user"] --> W["Django web application"]
    W --> DB[("Relational database")]
    W --> FS[("Static, media and research artifacts")]
    W --> R["Redis / Celery when enabled"]
    W --> DOCS["MkDocs and in-app documentation"]
```

The primary UI is server-rendered Django with incremental HTMX interactions. Development uses SQLite by default. PostgreSQL and a production-grade runtime are target M0-S concerns unless separately verified in a deployment.

## Current application containers

```mermaid
flowchart TB
    subgraph APP["Django project: mmportal"]
        CL["clinic\nstructured records and navigation"]
        SI["simulator\nmechanistic scenarios and design prototypes"]
        TW["twin_engine\nstate, lineage, runs and research diagnostics"]
        CH["chemtools\noptional molecular utilities"]
        DV["docs_viewer\nin-app documentation"]
        GV["governance\nmachine-readable release claims"]
    end

    CL --> TW
    TW --> SI
    SI --> CH
    GV --> CL
    GV --> SI
    GV --> TW
```

## Bounded responsibilities

### `clinic`

- structured patient, assessment and treatment entities;
- patient list/detail and clinical-history surfaces;
- ownership fields and selected API filtering;
- current product-authentication entry points;
- temporary governed smoke-user command.

Current limitation: object-level authorization is not yet consistently centralized across all web views and aggregate queries.

### `simulator`

- `MathematicalModel` ODE execution;
- scenario and `SimulationAttempt` workflows;
- PK/PD presets and time-resolved dose functions;
- model-relative KPI summaries;
- experimental therapy-design code.

Current limitation: core biological parameters and several endpoint definitions are heuristic research abstractions.

### `twin_engine`

- persistent `PatientTwinState` and parent/child lineage;
- computational-input identity;
- longitudinal labs, interruptions and adverse events;
- observation model and residuals;
- calibration infrastructure;
- mechanistic counterfactual runs;
- uncertainty, sensitivity, robustness and backtesting prototypes;
- immutable run manifests, artifact hashes, comparability and invalidation;
- Simple Research View, Research Cockpit and Developer Console.

Current limitation: infrastructure does not establish identifiable or externally validated patient-specific prediction.

### `chemtools`

- optional RDKit-dependent and molecular-analysis utilities;
- asynchronous jobs where configured;
- generated artifacts.

Current limitation: molecular utilities are exploratory and are not evidence of anti-myeloma efficacy.

## Current scientific run flow

```mermaid
sequenceDiagram
    autonumber
    actor U as Research user
    participant C as clinic/twin data
    participant T as twin_engine
    participant S as simulator
    participant A as artifact store

    U->>T: Select governed Twin and intervention
    T->>C: Validate patient access, input identity and readiness
    T->>S: Execute mechanistic baseline/alternative simulation
    S-->>T: Trajectory, KPIs and predicted biomarker projection
    T->>T: Add toxicity, uncertainty and comparison diagnostics
    T->>A: Write hashed trajectory/report artifacts
    T->>C: Persist immutable run manifest and artifact records
    T-->>U: Render research-only result and limitations
```

## Current data and identity flow

```mermaid
flowchart LR
    S["Source assertion / private evidence"] --> D["Versioned structured dataset"]
    D --> I["Idempotent import"]
    I --> CI["Content-addressed computational input"]
    CI --> TS["Lineage-bound Twin state"]
    TS --> RM["Immutable research run manifest"]
    RM --> AR["Hashed artifacts"]
```

Private source documents and direct identifiers stay outside Git. The repository stores contracts and de-identified/synthetic material only.

## Current authentication boundary

```text
public: documentation and selected simulator surfaces
protected: clinic, research, simulator-management and API surfaces
configured login target: Django admin login
normal product login/role contract: incomplete
```

The target least-privilege role architecture is tracked in M0-S issue `#8`; it is not current behavior.

## Current deployment and CI boundary

- deterministic dependency lock exists;
- several GitHub Actions workflows exist;
- image-build and deployment automation exist;
- branch-protected aggregate scientific CI is not yet verified;
- production security, PostgreSQL, storage, quotas, observability and recovery must be evidenced separately under M0-S.
