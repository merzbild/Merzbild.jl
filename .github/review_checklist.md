### Review checklist

This checklist is meant to assist creators of PRs (to let them know what reviewers will typically look for) and reviewers (to guide them in a structured review process). Items do not need to be checked explicitly for a PR to be eligible for merging.

#### AI use
- [ ] Any use of AI tools is disclosed.
- [ ] A self-assessment on [https://www.visidata.org/blog/2026/ai/](0-10 scale) is provided (see CONTRIBUTING.md for more details).

#### Code quality
- [ ] The code can be understood easily.
- [ ] Newly introduced names for variables etc. are self-descriptive and consistent with existing naming conventions.
- [ ] There are no redundancies that can be removed by simple modularization/refactoring.
- [ ] There are no leftover debug statements or commented code sections.
- [ ] The code adheres to the coding convention laid out in CONTRIBUTING.md.

#### Documentation
- [ ] New functions and types are documented with a docstring; this is link in the public or private API documentation.
- [ ] Relevant publications are referenced in docstrings.
- [ ] Inline comments are used to document longer or unusual code sections.
- [ ] Comments describe intent ("why?") and not just functionality ("what?").
- [ ] Changes are documented in CHANGELOG.md.

#### Testing
- [ ] The PR passes all tests.
- [ ] New or modified lines of code are covered by tests.

#### Performance
- [ ] There are no type instabilities or memory allocations in performance-critical parts.
- [ ] If the PR intent is to improve performance, before/after are posted in the PR.

#### Verification
- [ ] The correctness of the code was verified using appropriate tests.
- [ ] If new equations/methods are added, a convergence test has been run and the results
  are posted in the PR.