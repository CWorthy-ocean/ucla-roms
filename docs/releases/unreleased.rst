.. _unreleased:

Unreleased
----------

.. note::
    This release is currently in development

Breaking Changes
~~~~~~~~~~~~~~~~

- N/A

New Features
~~~~~~~~~~~~

- N/A

Bug Fixes
~~~~~~~~~


- Add explicit NETCDFF paths to `make` linking (`#313 <https://github.com/CWorthy-ocean/ucla-roms/pull/313>`_)
- Fix missing cdr tracer param usage (`#317 <https://github.com/CWorthy-ocean/ucla-roms/pull/317>`_)
- Fix SSS test using coarse dims, which isn't supported for restoring (`#317 <https://github.com/CWorthy-ocean/ucla-roms/pull/317>`_)

Improvements
~~~~~~~~~~~~

- N/A

Miscellaneous
~~~~~~~~~~~~~

- Fixes incorrect backtick usage in the release notes, and updates the PR template scraper to use backticks correctly (`#315 <https://github.com/CWorthy-ocean/ucla-roms/pull/315>`_)
- Adds github action to finalize release notes ahead of a new release (`#315 <https://github.com/CWorthy-ocean/ucla-roms/pull/315>`_)
- CI: Pins parts of the gnu compiler stack (`#316 <https://github.com/CWorthy-ocean/ucla-roms/pull/316>`_)
- CI: xfail on hash mismatches, hard fail on other errors (`#316 <https://github.com/CWorthy-ocean/ucla-roms/pull/316>`_)
- CI: update roms-tools to fix test construction failures (`#316 <https://github.com/CWorthy-ocean/ucla-roms/pull/316>`_)
- Standardize trailing `/` convention for Makedefs (`#313 <https://github.com/CWorthy-ocean/ucla-roms/pull/313>`_)
- Allow additional link flags to be passed via USER_LDFLAGS (`#313 <https://github.com/CWorthy-ocean/ucla-roms/pull/313>`_)
- Update release notes finalizer to remove sections with nothing in them (`#319 <https://github.com/CWorthy-ocean/ucla-roms/pull/319>`_)
- Update release notes updater to handle unbulleted content (`#319 <https://github.com/CWorthy-ocean/ucla-roms/pull/319>`_)
