"""Entry point of the own generator bank: ALL_TASKS is the list of Task (see common.py).

Example:
    from registry import ALL_TASKS, get
    task = get("C05_qft")
    qc = task.styles[0].fn(n=5)
"""
from tasks_classic import TASKS as _C
from tasks_qcoder import TASKS as _Q

ALL_TASKS = list(_C) + list(_Q)


def get(task_id):
    for t in ALL_TASKS:
        if t.id == task_id:
            return t
    raise KeyError(task_id)
