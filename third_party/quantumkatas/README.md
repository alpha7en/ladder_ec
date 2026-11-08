# Qiskit QuantumKatas material

Everything in this folder is derived from Qiskit QuantumKatas and is licensed under the Creative Commons
Attribution-NonCommercial-ShareAlike 4.0 International License (CC BY-NC-SA 4.0), see [LICENSE](LICENSE).
It is not covered by the MIT licence of the rest of the repository.

- Source: https://github.com/qiskit-community/qiskit-quantumkatas, commit
  `f93049650950bab00a8a1681f304a38fb0c02ad3` (2026-10-20).
- Copyright (c) 2018-2024 Microsoft Corporation; copyright (c) 2025-2026 contributors to the Qiskit-QuantumKatas
  project. Based on Microsoft's QuantumKatas (https://github.com/microsoft/QuantumKatas).
- Paper: J. Cruz-Benito, I. Faro, "Qiskit QuantumKatas: Adapting Microsoft's Quantum Computing Exercises for LLM
  Evaluation", arXiv:2605.27210, 2026.

| File | What it is | Changes |
|---|---|---|
| `dataset.jsonl` | `dataset/qiskit_quantumkatas.jsonl` of the source (350 tasks) | none |
| `wrappers/<task>__A.py` | reference solution of a kata plus a wrapper `solve(n)` or `solve(n, a)` | wrapper appended |
| `wrappers/<task>__B.py` | the same without the final `return qc` of each function | final returns removed, wrapper appended |
| `llm_solutions.jsonl` | solutions written by a language model for 52 katas; each starts with the task statement of the kata as a comment | new code under the quoted statement |

The wrappers are written by `experiments/external/quantumkatas/run_katas.py`; the solutions are evaluated by
`experiments/external/quantumkatas/llm_katas.py` (Sec. 7.5 of the paper, Table 10).
