# Release-note fragments

A pull request that changes something users could notice adds a release note: a file in this directory named `+<name>.<type>.rst`, where `<name>` is any short label, such as `+river-volume.bugfix.rst`. Commit it with your code. The `+` means the pull request number is filled in at release time, from the merge commit, so the note links to your pull request. If you already know the number you can use it instead: `<PR number>.<type>.rst`. At release time the notes are collected into `docs/releases/<version>.rst` and the files here are deleted, so don't edit `docs/releases/` yourself.

The type decides the section the note appears under:

| Type          | Section           | For |
|---------------|-------------------|-----|
| `breaking`    | Breaking Changes  | a setup that worked now fails or needs an edit: a renamed or removed namelist key or CPP key, an input file that is now rejected, an output file that code reading it (such as C-Star) can no longer use |
| `output`      | Changes to Output | the same setup still builds and runs, but its results or the files it writes differ: values, variables, or the number of files |
| `deprecation` | Deprecations      | something that still works but will be removed |
| `feature`     | New Features      | new capabilities |
| `bugfix`      | Bug Fixes         | code that did not do what it should, and now does |
| `improvement` | Improvements      | better behavior, messages or performance of existing code |
| `misc`        | Miscellaneous     | documentation, CI, build and test changes worth mentioning |

Write the note for someone reading the release notes, not for the reviewer: what changed for them and anything they need to do, in one line. Details for the reviewer go in the pull request description. If the reader needs more, add `- ` sub-bullets on the following lines. The notes are reStructuredText, so code, namelist keys and file names go in double backticks; outside them a name that ends in an underscore, such as `ub_`, is read as a link and fails the check on the pull request. For example, `changes/+avg-temp-salt.feature.rst`, added by pull request 413:

```rst
``wrt_avg_temp`` and ``wrt_avg_salt`` in ``TS_OUTPUT_SETTINGS`` add ``temp`` and ``salt`` to the averages file (with ``wrt_file_avg``).
- Both are off by default, and a namelist without them runs as before.
```

becomes, in the release notes:

```rst
New Features
~~~~~~~~~~~~

- ``wrt_avg_temp`` and ``wrt_avg_salt`` in ``TS_OUTPUT_SETTINGS`` add ``temp`` and ``salt`` to the averages file (with ``wrt_file_avg``). (`#413 <https://github.com/CWorthy-ocean/ucla-roms/pull/413>`_)

  - Both are off by default, and a namelist without them runs as before.
```

Each note is its own file, so a pull request with several notes adds several files, such as `+inflow-flux.output.rst`, `+inflow-flux.bugfix.rst` and `+obc-docs.misc.rst`. Two notes of the same type need different names. A change that both breaks setups and changes results gets one note of each type. A deprecation note reads ``Deprecated: <what>, removed in <release>; use <replacement>.``

A pull request with nothing for users to read about, such as a test-only change or a typo fix, needs no note; ask a maintainer to add the `no-changelog` label. The `towncrier check` job on the pull request fails until it has a note or the label, and it also fails on a note that is not valid reStructuredText.

To preview the notes as they will appear, install [towncrier](https://towncrier.readthedocs.io) (`mamba install -c conda-forge towncrier`) and run this from the repository root:

```bash
towncrier build --draft --version Unreleased
```

Notes named with a `+` show up without their pull request link until the release. `ci/README` describes how a release is made.
