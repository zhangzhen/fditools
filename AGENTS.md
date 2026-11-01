# fditools

Rust CLI `fdi` for sequencing-data processing (FASTQ, BAM), built on the crate's library
`fditools_lib`.

## Commands

- Build: `cargo build --release` (binary at `target/release/fdi`)
- Test: `cargo test`
- Type-check: `cargo check --all-targets`

## Docs

- [Agent setup](.agents/docs/agent-setup.md): read when adding a doc, skill or rule.

## Planning

Save plans, requirements and investigations in `.agents/plans/<slug>/<slug>.<TYPE>.md`
(gitignored; e.g. `some-task.PLAN.md`), prefixing the slug with the ticket id when there is one.
