# Code Review: LLVM IR Optimization and SIMD Vectorization

## Command: `/codeReviewLLVM`

When this command is invoked, perform a comprehensive review of Swift code for LLVM-friendly patterns and SIMD optimization opportunities, with a focus on Linux/x86_64 targets.

## Review Checklist

### 1. SIMD Storage Types
- [ ] Check if using `SIMD2<T>`, `SIMD3<T>`, or `SIMD4<T>` for storage
- [ ] Verify direct SIMD operations instead of scalar operations
- [ ] Ensure storage is marked `@usableFromInline` when used in `@inlinable` functions

**Good:**
```swift
@usableFromInline
internal var storage: SIMD2<T>

let result = lhs.storage + rhs.storage
```

**Bad:**
```swift
var real: T
var imaginary: T

let result = Complex(real: lhs.real + rhs.real, imaginary: lhs.imaginary + rhs.imaginary)
```

### 2. Inlining Annotations
- [ ] All hot-path functions marked `@inlinable`
- [ ] All computed properties on critical path marked `@inlinable`
- [ ] Initializers that should inline marked `@inlinable`

**Check for:**
- Mathematical operations (+, -, *, /)
- Normalization functions
- Magnitude calculations
- Type conversions

### 3. SIMD Intrinsics Usage
- [ ] Use `simd_length()` instead of manual sqrt of dot product
- [ ] Use `simd_length_squared()` instead of manual dot product
- [ ] Use `simd_normalize()` instead of manual normalization
- [ ] Use `__sincos()` / `__sincosf()` for simultaneous sin/cos

**Good:**
```swift
let magnitude = simd_length(storage)
let normalized = simd_normalize(storage)
var s: T = 0, c: T = 0
__sincos(angle, &s, &c)
```

**Bad:**
```swift
let magnitude = sqrt(x*x + y*y)
let normalized = self / magnitude
let s = sin(angle)
let c = cos(angle)
```

### 4. Memory Layout
- [ ] Structs have appropriate alignment for SIMD
- [ ] No unnecessary padding
- [ ] Storage fields are contiguous
- [ ] Public fields use computed properties backed by SIMD storage

### 5. Control Flow
- [ ] Minimize branches in hot paths
- [ ] Use select/conditional moves instead of if/else when possible
- [ ] Avoid complex control flow in tight loops

**Good:**
```swift
// Branchless
let result = condition ? trueValue : falseValue
```

**Bad:**
```swift
var result: T
if condition {
    result = expensiveComputation1()
} else {
    result = expensiveComputation2()
}
```

### 6. Fast Math Opportunities
- [ ] Check if operations can use fast math flags
- [ ] Ensure no precision-critical code prevents optimization
- [ ] Look for reassociation opportunities

### 7. Loop Vectorization
- [ ] Loops with simple, regular access patterns
- [ ] No loop-carried dependencies
- [ ] Trip count known or easily computed
- [ ] Operations suitable for vectorization

### 8. Function Call Overhead
- [ ] No unnecessary function calls in tight loops
- [ ] Generic functions are specialized where possible
- [ ] No dynamic dispatch on hot paths

## Analysis Steps

When reviewing code:

1. **Identify Hot Paths**: Locate performance-critical operations
   - Mathematical operations
   - Frequently called functions
   - Inner loops

2. **Check SIMD Opportunities**:
   ```
   Look for patterns like:
   - x*x + y*y → simd_length_squared(SIMD2(x,y))
   - sqrt(x*x + y*y) → simd_length(SIMD2(x,y))
   - Manual vector operations → SIMD operations
   ```

3. **Verify Inlining**:
   - All operations <10 lines should be `@inlinable`
   - Critical path must have complete inline chain

4. **Recommend Improvements**:
   - Specific code changes
   - Expected performance impact
   - LLVM IR patterns to verify

## Output Format

Provide review in this structure:

```markdown
## LLVM/SIMD Code Review

### ✅ Strengths
- List well-optimized patterns found
- Note good use of SIMD intrinsics
- Highlight proper inlining

### ⚠️  Optimization Opportunities
For each issue:
1. **Location**: File:Line
2. **Issue**: What's suboptimal
3. **Current Code**: Show current implementation
4. **Recommended**: Show optimized version
5. **Impact**: Expected performance improvement
6. **LLVM IR Check**: How to verify in IR

### 🔍 Verification Commands
Provide specific commands to verify optimizations:
- LLVM IR generation
- Vector instruction grep patterns
- Assembly inspection

### 📊 Expected SIMD Score
- Current estimated score: X/10
- After improvements: Y/10
```

## Example Review Output

```markdown
## LLVM/SIMD Code Review for Complex.swift

### ✅ Strengths
- Using SIMD2<T> storage ✓
- All arithmetic operators marked @inlinable ✓
- Using simd_length() for magnitude ✓

### ⚠️ Optimization Opportunities

#### 1. Phasor initialization can use @inlinable
**Location**: Complex.swift:96
**Issue**: Phasor init not marked @inlinable
**Impact**: Prevents cross-module optimization
**Fix**: Add @inlinable attribute

#### 2. Conjugate using scalar operations
**Location**: Complex.swift:47
**Issue**: Manual negation instead of SIMD
**Current**:
```swift
return Complex(real: storage.x, imaginary: -storage.y)
```
**Recommended**:
```swift
return Complex(vector: storage * SIMD2(1, -1), isNormalized: _isNormalized)
```
**Impact**: Better vectorization, fewer operations
**LLVM IR Check**: Look for `fmul <2 x T>` instead of separate ops

### 🔍 Verification Commands
```bash
# Generate IR
swiftc -emit-ir -O -target x86_64-unknown-linux-gnu Complex.swift -o Complex.ll

# Check for vector operations
grep '<2 x ' Complex.ll | wc -l

# Look for SIMD intrinsics
grep 'llvm.fmuladd\|llvm.sqrt' Complex.ll
```

### 📊 SIMD Optimization Score
- Current: 8/10 (Excellent base, minor improvements possible)
- After fixes: 9/10
```

## Focus Areas for Linux/x86_64

- Prioritize AVX/AVX2 vectorization (256-bit)
- Check for SSE fallbacks
- Verify alignment for cache line efficiency
- Look for opportunities to use FMA (fused multiply-add)
- Ensure no macOS-specific patterns (vDSP, Accelerate)

## Red Flags

Immediately flag these patterns:
- ❌ Manual sqrt of sum of squares (use `simd_length`)
- ❌ Missing `@inlinable` on <10 line functions
- ❌ Scalar storage for vector quantities
- ❌ Separate sin/cos calls (use `__sincos`)
- ❌ Complex branches in arithmetic operations
- ❌ Non-SIMD type conversions in hot paths
