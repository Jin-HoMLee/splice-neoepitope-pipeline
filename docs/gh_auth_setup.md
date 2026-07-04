# Local GitHub CLI Authentication Setup

Best-practice setup for `gh` authentication on a local development machine.

The goal is to stop relying on a long-lived classic PAT hardcoded in a plaintext dotfile, and instead make the system keyring (macOS Keychain) the primary interactive credential.

This document is the Developer-side advisory artifact for [Issue #971](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/971).
The account-level actions it points to (the browser OAuth flow, editing `~/.zshenv`, revoking tokens) are **human-only** and are listed as a checklist at the end.

## Why change the current setup

The prior setup is the anti-pattern on three counts, each confirmed against GitHub's own guidance (sourced below):

1. **Long-lived classic PAT with broad `repo` scope.**
   GitHub recommends fine-grained, minimally-scoped, expiring tokens.
   A leaked broad classic token can touch every repo the account owns.
2. **Plaintext in a dotfile** (`export GH_TOKEN=...` in `~/.zshenv`).
   Readable by any process running as the user; the Keychain avoids on-disk plaintext.
3. **`GH_TOKEN` always-on for interactive use.**
   The `GH_TOKEN` env var shadows the Keychain on every `gh` call and is intended for headless automation, not day-to-day interactive work.

This is not theoretical: the token expired mid-session on 2026-07-03 and broke all `gh` operations with 401s partway through a board restructure, and two token values appeared in that session's command output.

## Target setup

Two credential paths, used for different purposes:

| Use | Credential | Storage |
|-----|-----------|---------|
| **Interactive** (day-to-day `gh` in a terminal) | OAuth token from `gh auth login` web flow | macOS Keychain (managed by `gh`) |
| **Automation only** (headless scripts, CI-like local runs) | fine-grained, expiring PAT | injected into `GH_TOKEN` at runtime from the Keychain, never a plaintext dotfile |

Key point: do **not** keep a plaintext `GH_TOKEN` in `~/.zshenv` for interactive use.
When `GH_TOKEN` is set, it overrides the Keychain credential for every `gh` invocation, which is exactly the shadowing that caused the mid-session breakage.

### Interactive: Keychain OAuth via `gh auth login`

```bash
# Remove the always-on env token first (see the human checklist), then:
gh auth login --hostname github.com --git-protocol https --web
```

The web flow stores the OAuth token in the macOS Keychain.
Its default scopes do not include `project`, and our board is a user Projects v2 board (user project 9).
Add `project` (which alone grants Projects v2 read/write on user- or org-owned boards) plus `read:org` (a `gh auth login` default that AC-4 mandates; it covers org membership, not the board ops themselves):

```bash
gh auth refresh --hostname github.com --scopes project,read:org
```

Verify (should show `project`, `read:org`, `repo`, and no missing-scope warning):

```bash
gh auth status
```

### Automation: fine-grained expiring PAT (only if genuinely needed)

Create the token at <https://github.com/settings/personal-access-tokens> with:

- **Resource owner:** `Jin-HoMLee`
- **Repository access:** only the repos automation touches (at minimum `splice-neoepitope-pipeline`; add the personas repo if a routine drafts there)
- **Expiration:** set an expiry (e.g. 90 days); a non-expiring token is the anti-pattern being removed
- **Repository permissions:**
  - Contents: Read and write
  - Issues: Read and write
  - Pull requests: Read and write
  - Metadata: Read-only (mandatory, auto-selected)
  - Workflows: Read and write (only if a routine edits `.github/workflows/`)
- **Account permissions:**
  - Projects: Read and write (this is what drives board operations on the user project; the classic `project` scope equivalent for fine-grained tokens)

Note: GitHub's own `gh` guidance favors supplying a fine-grained PAT via `GH_TOKEN` rather than storing it through `gh auth login`, because a fine-grained token's resource scoping can produce confusing behavior when the Keychain path assumes a broadly-scoped credential.
So the fine-grained token belongs in `GH_TOKEN`, but injected at runtime, not written to a dotfile.

### Keychain-lookup snippet (for a retained automation token)

Store the fine-grained token once in the Keychain:

```bash
security add-generic-password -U -a "$USER" -s gh-automation-token -w
# (paste the token at the prompt; it is not echoed and not stored on disk in plaintext)
# -U updates the item in place if it already exists, so this same command is safe to re-run on rotation
```

Then, in an automation shell only, inject it at runtime instead of exporting a literal:

```bash
# put in the automation script, NOT in ~/.zshenv:
export GH_TOKEN="$(security find-generic-password -a "$USER" -s gh-automation-token -w)"
```

Rotate by re-running the same `security add-generic-password -U ...` command with the new token after the old one is revoked (the `-U` flag updates the existing item in place); nothing in version control or dotfiles needs to change.

## Verification

After the interactive setup, both of these should succeed with the Keychain credential (no `GH_TOKEN` in the environment):

```bash
gh auth status                      # shows project, read:org, repo; no missing-scope warning
scripts/board_open_items.py --role developer   # a paginated board query that exercises the project scope
```

## Human-only checklist (Jin-Ho)

An agent cannot perform these; they are account-level or local-machine actions.

- [ ] **Rotate the two exposed classic PATs first (time-sensitive).**
  Revoke both at <https://github.com/settings/tokens>: the one that already expired on 2026-07-03, and the classic `ghp_` token still active as the current `gh` credential. Identify them by their name/note and creation date on that page; no token value is reproduced here.
- [ ] Run the `gh auth login --web` browser flow to store an OAuth token in the Keychain.
- [ ] Run `gh auth refresh --scopes project,read:org` to add the board scopes.
- [ ] Remove the `export GH_TOKEN=...` line from `~/.zshenv`.
- [ ] If automation genuinely needs a token, create the fine-grained PAT per the spec above and store it in the Keychain (not a dotfile).
- [ ] Confirm `gh auth status` shows `project`, `read:org`, `repo` with no missing-scope warning.

## Sources

- [Managing your personal access tokens](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens)
- [Introducing fine-grained personal access tokens (GitHub Blog)](https://github.blog/security/application-security/introducing-fine-grained-personal-access-tokens-for-github/)
- [`gh` environment manual](https://cli.github.com/manual/gh_help_environment)
- [`gh auth login` manual](https://cli.github.com/manual/gh_auth_login)
- [cli/cli #9461 - minimum required v2 fine-grained PAT scopes](https://github.com/cli/cli/issues/9461)
