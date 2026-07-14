---
chunk: 07-packaging-and-linux
status: pending
depends_on: []
audit: ../review-for-fixes/2026-07-11-swift-audit.md §F8, §F10, §F11; decision: Linux is a required target
last_updated: 2026-07-11
semver: 0.0.1
author: Nicholas Bergantz
---

# 07 — Packaging: pins, platforms, Linux target

**Deliverable:** reproducible dependency pinning, matched platform list,
and a Linux verification target.

## Files

- Edit: `Package.swift`, `makefile`
- Possibly edit: whatever the Linux build flushes out (small
  `#if canImport` patches only; anything larger → stop and report)

## Design constraints

1. **F8 (pin):** `branch: "dev"` on spmFoundationTools becomes an
   `exact:`/`from:` pin of a tagged release. **Human prerequisite:** a tag
   must exist on spmFoundationTools (its Makefile has `make tag`) —
   this chunk BLOCKS until the tag is cut; note the required tag in the
   chunk status if blocked. Interim dev workflow: document the local
   path-override pattern (`swift package edit` or a local `.package(path:)`
   swap) in the readme's development section instead of the branch pin.
2. **F10:** remove `depermaid` from dependencies (verify `make mermaid`
   behavior first, same rule as the foundation plan's chunk 06).
3. **F11 (platforms):** widen `platforms:` to match spmFoundationTools
   (`macOS .v14, iOS .v16, tvOS .v16, watchOS .v9`) — nothing in this
   target is macOS-specific; the compiler will confirm.
4. **Linux:** add the same containerized target as the foundation plan:

```make
linux-test:  ## Build & test in a Linux container (requires Docker)
	docker run --rm -v "$(PWD)":/pkg -w /pkg swift:6.1 \
		swift test --package-path .
```

   Note: math has no Darwin-only APIs (`__sincos` lives in foundation), so
   this should pass once foundation's chunks 01/02 are released in the
   pinned tag — sequence accordingly.

## TDD steps

1. `make build test` before/after each Package.swift edit;
   `make linux-test` runs and its outcome is recorded (pass, or blocked on
   the foundation tag — say which).

## Acceptance criteria

- [ ] `grep -n "branch:" Package.swift` → no hits
- [ ] `grep -n "depermaid" Package.swift` → no hits
- [ ] Platform list matches foundation's
- [ ] `make build test` passes; `make linux-test` outcome recorded

## Out of scope

Foundation-repo changes (its own plan); CI wiring.
