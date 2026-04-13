# Copilot Agent Instructions – R Package Development

## 1. Scope and Role

You are assisting with the development of **research‑oriented R packages**, many of
which integrate **Rcpp / C++** components for performance. Your primary goal is to
produce **correct, minimal, CRAN‑compliant R code** that supports reproducible
scientific workflows.

You are **not** writing exploratory scripts unless explicitly asked.

---

## 2. Hard Constraints (Non‑Negotiable)

- All changes **must pass `R CMD check --as-cran`** without warnings or notes.
- **Never** introduce:
  - global side effects or modification of the user’s environment
  - file I/O inside package functions (unless the function’s explicit purpose is
    reading/writing a specified file)
  - assumptions about interactive sessions
  - system‑specific paths (e.g., `C:\...` or `/home/user/...`)
- Linux‑only commands are acceptable; Windows workarounds are **not required**.

If a request would violate any of these constraints, **explain why before proceeding**.

---

## 3. R Package Structure Rules

Respect the standard R package layout:

- `R/`
  - Contains **only** function definitions.
  - **No executable code** at the top level (i.e., no `library()` calls, no `options()`).
- `man/`
  - Generated **exclusively** by **roxygen2**.
  - Never edit `.Rd` files manually unless explicitly instructed.
- `src/`
  - C++ source files (`*.cpp`, `*.h` / `*.hpp`).
  - Rcpp attributes (`// [[Rcpp::export]]`, `// [[Rcpp::depends(...)]]`) are required.
- `tests/`
  - Use **testthat** for all unit tests.
- `inst/`
  - Only for shipped data, templates, or examples.
- `vignettes/`
  - Long‑form documentation only.

Do **not** place:
- analysis scripts in `R/`
- helper scripts outside the proper package structure
- data‑loading code inside functions (use lazy‑data or pass as arguments)

---

## 4. Coding Standards (R)

- **Prefer base R** unless a dependency is **explicitly justified** (e.g., `checkmate`
  for validation, `Rcpp` for performance).
- **Always use explicit namespaces**: `stats::optim()`, `utils::head()`, `checkmate::assert_numeric()`.
- **Input validation** must be performed with **checkmate**.
  - Use `checkmate::assert_*()` early in every exported function.
  - Internal helpers may rely on the caller’s validation.
- **Avoid non‑standard evaluation** unless the function is specifically designed for
  interactive use (e.g., `dplyr` verbs); otherwise use standard evaluation.
- **Prefer clear control flow** over clever but obscure vectorisation.
- Every exported function must:
  - validate all inputs
  - fail early with an informative error message
  - be **deterministic** (no random behaviour without `set.seed()`)

Internal helper functions should **not** be exported unless required for advanced
users. Use `@keywords internal` in roxygen2.

### 4.1. C++ / Rcpp Integration

- Rcpp code must be **header‑only** where possible (e.g., `featured_space.hpp`).
- Use **`Rcpp::XPtr`** to pass compiled C++ objects back to R (as seen in
  `maxentcpp`, `nicher`, `ucminfcpp`).
- Prefer **`RcppParallel`** for embarrassingly parallel loops; always reset the
  thread count after use.
- Finite‑difference gradients in C++ should respect user‑supplied `gradstep`
  vectors (see `ucminfcpp` and `niche_obj.cpp`).

### 4.2. Optimization Workflow

- For general unconstrained optimisation, use **`ucminf::ucminf`** (or
  `ucminfcpp::ucminf_xptr` when the objective is compiled C++).
- Provide a `control` argument that merges user values with sensible defaults.
  - Example: `control = list(grad = "central", gradstep = c(1e-6, 1e-8), maxeval = 500)`.
- Use **`do.call(ucminf_control, control)`** to build the final control list.

---

## 5. Documentation Rules (roxygen2)

- Every exported function **must** have:
  - `@param` entries for **all** arguments
  - `@return` description (even if `NULL`)
  - at least **one runnable example**
- Examples must:
  - run **quickly** (typically < 5 seconds)
  - use `set.seed()` when randomness is involved
  - avoid file I/O and large data downloads
- **Internal functions** must be marked with `@keywords internal`.
- For packages that use compiled code, include:
  - `@useDynLib <pkgname>, .registration = TRUE`
  - `@importFrom Rcpp sourceCpp`
- Use **`@seealso`** to cross‑reference related functions.

---

## 6. Testing & Reproducibility

- **Testthat** is the only accepted testing framework.
- Prefer **small, deterministic unit tests**.
- Always use `set.seed()` when generating random data.
- Explicitly state numerical tolerances and justify them.
- Tests should cover:
  - edge cases (empty inputs, `NA`, `Inf`, extreme values)
  - failure modes (informative error messages)
  - invariants (e.g., weights sum to 1, matrix orthogonality)
- For C++ functions, include **cross‑language consistency tests** comparing
  pure‑R and C++ implementations (see `test-like_neg_ltsgr_cpp_vs_r.R`) in
  preliminary versions and after the test pass remove the older r code and
  keep only the c++ version.

---

## 7. Dependency Policy

- **New dependencies require explicit justification.**
- Prefer **Suggests** over **Imports** for optional functionality (e.g., `future`
  for parallel optimisation, `pomp` for Sobol’ sequences).
- Do **not** add a dependency solely for a convenience function if a base‑R
  alternative exists.
- When proposing a new dependency, always state:
  - why it is needed
  - whether a base‑R alternative is viable
- Avoid `tidyverse` dependencies unless the package is specifically designed for
  interactive data manipulation (your packages are **not**).

---

## 8. Performance and Numerical Stability

- **Avoid unnecessary allocations** inside loops; pre‑allocate vectors/matrices.
- Prefer **stable algorithms** over fast but fragile ones (e.g., log‑sum‑exp for
  likelihoods).
- **Explicitly guard against `NaN` / `Inf`** by returning a large finite penalty
  (`1e300`) so the optimiser can recover.
- Use `constexpr` for constants in C++ (e.g., `OPTIM_PENALTY`, `MIN_KDE_WEIGHT`).
- For C++ code, favour **Eigen** for linear algebra and ensure row‑major / column‑major
  conventions are clearly documented.

---

## 9. Interaction Rules

- **Do not ask clarifying questions unless the task is truly ambiguous.**  
  If multiple designs are possible:
  - present the **recommended approach first**
  - briefly mention one alternative with a one‑sentence trade‑off
- **Prefer minimal, focused answers** – no bonus sections or tangents.
- When uncertain about a design choice, **state the uncertainty clearly** and
  provide a conditional path.

---

## 10. Output Style

- Use **clean Markdown**.
- **Never** place code inside Markdown tables.
- Use **fenced code blocks** with the correct language tag (````r` or ````cpp`).
- Avoid verbosity unless explicitly requested.
