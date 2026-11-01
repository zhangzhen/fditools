# Agent setup

How this repo is wired for coding agents, and how to hook up yours.

## Layout

- `AGENTS.md`: the canonical entry file and docs index. Keep it lean; everything else goes into a
  doc linked from it.
- `.agents/docs/`: agent-facing docs.
- `.agents/skills/<name>/`: project skills, each with a committed relative symlink
  `.claude/skills/<name>` → `../../.agents/skills/<name>`.
- `.agents/plans/`: per-developer plans, gitignored.

## Hooking up an agent

- **Claude Code** v2.1.277+ reads `AGENTS.md` directly; the repo has no `CLAUDE.md` on purpose.
  Any `CLAUDE.md`, `.claude/CLAUDE.md` or `CLAUDE.local.md` here or in a parent directory makes
  Claude read it **instead**. That includes your gitignored `CLAUDE.local.md`: open it with
  `@AGENTS.md`, or set **Project instructions** to `claude-md-and-agents-md`.
- Sessions that can't load `AGENTS.md` (Claude Code < 2.1.281 on Bedrock, Vertex AI, a gateway or
  with telemetry off; the `agents-md` plugin disabled; **Project instructions** set to
  `claude-md`) need a root `CLAUDE.md` holding only `@AGENTS.md`.
- **Any other agent**: point its entry-file mechanism at `AGENTS.md` (an include where supported,
  else a copy kept in sync), or re-run `/agentify-project` to research and wire it.

## Adding docs, skills and rules

- Doc: put it in `.agents/docs/` and link it from `AGENTS.md` with a plain markdown link and a
  one-line "read when…" hook. Never `@`-import it: imports load at launch and defeat progressive
  disclosure.
- Skill: put it in `.agents/skills/<name>/` and add the symlink above.
- Hook test: after moving content out of `AGENTS.md`, give a fresh session a task needing that doc
  and check it opens it. A miss means a weak hook; sharpen it (e.g. via `/self-improve`).

## Windows

Symlinks need `git config core.symlinks true` plus Developer Mode (or admin rights); otherwise
they check out as plain text files.
