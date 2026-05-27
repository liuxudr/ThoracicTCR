"""Per-stage Typer sub-apps for the thoracictcr CLI.

Each module exposes a single Typer `app` which is mounted by ``thoracictcr.cli``.
The CLI is the user-facing entry point (``thoracictcr``); the ``scripts/`` shell
and Python files remain as thin wrappers that delegate to these sub-apps so the
pipeline is reproducible from either side.
"""
