# spmMathTools — Python Backport Considerations

When backporting to Python, key architectural differences to plan for:

- Swift's value types (structs) → Python classes (all reference). Need explicit `copy.deepcopy()` or `dataclasses` with `frozen=True`
- Swift operator overloads → Python `__add__`, `__mul__`, etc. (direct mapping)
- Swift generics (`<T: BinaryFloatingPoint>`) → Python generics (`Generic[T]` with `TypeVar`) or `numpy` dtype parameterization
- `Sendable` concurrency safety → no direct Python equivalent; consider `multiprocessing`-safe patterns
- Polynomial class hierarchy → Python can use the same pattern or simplify with `numpy.poly1d` / `numpy.polynomial`
- OTG module → direct port feasible; Python Ruckig bindings already exist in `debugSupport/` for validation
