# Third-party material

Licences were checked on 2026-11-08 against the upstream repositories and the wheel at the versions below.

| Material | Version | In this repository | Licence |
|---|---|---|---|
| Qiskit QuantumKatas (qiskit-community/qiskit-quantumkatas) | commit `f93049650950bab00a8a1681f304a38fb0c02ad3` | `third_party/quantumkatas/`: the dataset, the wrapped reference solutions, language-model solutions that quote the task statements | CC BY-NC-SA 4.0 (`third_party/quantumkatas/LICENSE`); attribution in `third_party/quantumkatas/README.md` |
| QuanBench (GuoXiaoYu1125/Quanbench) | commit `562b6ee545dd32823c84cb173bd88e710533f2dc` | not included; fetched by `experiments/external/quanbench/fetch.sh`; per-function results in `data/runs/quanbench_translate.jsonl` | MIT (Copyright (c) 2026 QuanBench contributors) |
| QuanBench+ (JawadKotaichh/quanbench-plus) | commit `2dfd1a863b13d3762a734ed96742adb39e65e34b` | not included; fetched by the same script | MIT (Copyright (c) 2026 Jawad Kotaich) |
| MQT Bench (`mqt.bench` on PyPI) | 2.3.0 | not included; fetched by `experiments/external/mqtbench/fetch.sh`; results in `data/runs/mqtbench_*`, `data/external/mqtbench/classification.csv` | MIT (Copyright (c) 2023 - 2026 Chair for Design Automation, TUM) |
| AutoQ-Para (VeriFIT/AutoQ-Para) | commit `e3a4718166069086a4923b48dfda7ae2f944a72e` | not included; fetched and built by `experiments/autoq_para/fetch_and_build.sh`; only our build patch (`build_fixes.patch`) and our models (`autoq_models.cpp`) are here | MIT (Copyright (c) 2025 Michal Hečko) |
| Qiskit, NumPy, mqt.qcec | 2.5.2, 2.5.3, 3.11.0 | dependencies, installed by pip | Apache 2.0, BSD 3-Clause, MIT |

Our files `build_fixes.patch` and `autoq_models.cpp` modify and extend AutoQ-Para and are distributed under its MIT licence as well as ours.
