# Machines of the shipped runs

All runs in `data/runs/` were made in cloud containers of one type: 4 cores, 15 GB of memory, Python 3.13.18,
Qiskit 2.5.2, NumPy 2.5.3, mqt.qcec 3.11.0. Each run is one process on one core unless stated otherwise.
Times depend on the processor: the containers reported 2.10 GHz or 2.30 GHz Intel Xeon processors.

| Runs | Processor | Files |
|---|---|---|
| method, three variants, with retries | Intel Xeon @ 2.10 GHz | `method_main.jsonl`, `method_second.jsonl`, `method_no_quick_search.jsonl` |
| QCEC enumeration, `ncut` (start) and AutoQ-Para | Intel Xeon @ 2.10 GHz | `qcec_ncut.jsonl`, `autoq_para/` |
| QCEC enumeration (resumed), `ncut_min`, `witness` (start) | Intel Xeon @ 2.30 GHz | `qcec_ncut.jsonl`, `qcec_ncut_min.jsonl`, `qcec_witness.jsonl` |
| QCEC `witness` (resumed), instance checks | Intel Xeon @ 2.30 GHz | `qcec_witness.jsonl`, `instance_checks.jsonl` |
| step 4, language-model katas | Intel Xeon @ 2.10 GHz | `step4_ours.jsonl`, `step4_qcec.jsonl`, `llm_katas.jsonl`, `llm_katas_rerun.jsonl` |
| exact plan | Intel Xeon @ 2.10 GHz | `exact_plan.jsonl` |

The files without times (`katas_*.jsonl`, `quanbench_translate.jsonl`, `mqtbench_*`, `ladder_dims.csv`,
`own_generators_refusals.csv`) do not depend on the machine.
