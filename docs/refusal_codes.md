# Refusal codes

The translator refuses a generator outside the supported Python fragment, or a program outside the class of ladder
programs (Definition 3), with one of these codes. Full rules: Online Resource 1, Sect. S2, Table S2.

Kind: "class" means the code is in the fragment but the program is not a ladder program; "fragment" means the code
is outside the supported Python; "both" depends on the case.

| Code | Kind | Example |
|---|---|---|
| `branch_on_size` | class | `if n > N0`, `if n % 3 == 0` |
| `branch_on_loop_index` | class | `if i == N0`, `if i % 3 == 0` |
| `branch_on_param` | class | `if len(a) > 3` |
| `conditional_expression` | both | `range(n - 1 if n % 2 else n - 2)` |
| `nonaffine_expression`, `nonaffine_index` | both | `(i + 1) % n`, `range(n // 2)`, `2 ** j` |
| `nested_loop` | class | the double loop of QFT |
| `while_loop`, `recursion` | both | doubling trees, recursive QFT |
| `growing_arity` | class | `mcx`, `mcp` on a list whose length grows with n or i |
| `library_gate` | fragment | `compose`, `append`, `.inverse()` |
| `nonconst_angle` | class | `2 * acos(sqrt(1 / (n - i)))`, `pi / 2 ** i` |
| `two_registers`, `fold_index`, `stride_index` | class | indices `n + i`, `n - 1 - i`, `2 * i + 1` |
| `index_out_of_range`, `qubit_collision` | class | `cx(i, i + 1)` over `range(n)` on n qubits; `cx(0, i)` at i = 0 |
| `bit_index_negative`, `bit_index_out_of_range` | class | `a[i - 1]` at i = 0; `(a >> (i + 1)) & 1` at i = n - 1 |
| `bit_on_end_fixed` | class | a condition on the bit `a[n + c]` in a ladder body |
| `bit_const_in_down_ladder` | class | `for i in range(n - 1, 0, -1): if a[0]: qc.x(i)` |
| `fixed_overlaps_window`, `fixed_inside_ladder` | class | the window position `i - 1` meets the fixed qubit 0 at i = 1; `for i in range(n): if a[0]: qc.x(i)` |
| `alias_reassigned`, `alias_under_condition`, `loop_var_after_loop` | fragment | `t = t + 1`; assignment under `if a[i]`; loop counter used after the loop |
| `bit_param_misuse` | both | `range(a)`, `a[i]` outside a condition |
| `range_not_ladder` | both | `reversed(range(0, n, 2))`, a loop over an explicit list |
| `body_without_window` | class | `for _ in range(n): qc.x(0)` |
| `no_size_dependence`, `no_ladder` | both, fragment | a circuit of fixed size |
| `unknown_call`, `unknown_name`, `unknown_gate`, `unsupported_statement`, `unsupported_expression`, `extra_params`, `shadowed_builtin` | fragment | `getattr`, `try`, a default parameter, `QuantumRegister`, a redefined `range` |
| `early_return`, `nonint_constant`, `broadcast_multi`, `not_const`, `no_size_param` | fragment | `return` before the end; a float where an index is expected; no size parameter |
