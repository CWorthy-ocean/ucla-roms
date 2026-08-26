.. _unreleased:

Unreleased
----------

.. note::
    This release is currently in development

Breaking Changes
~~~~~~~~~~~~~~~~


- Existing child_bry files written with second-scale bry_time remain incompatible until regenerated or converted by dividing by 86400. (`#343 <https://github.com/CWorthy-ocean/ucla-roms/pull/343>`_)

New Features
~~~~~~~~~~~~

- N/A

Bug Fixes
~~~~~~~~~


- extract_data.F90 wrote bry_time from the seconds scalar time while labeling units as days, so nest OBC forcing failed to find the correct records / variables at runtime. (`#343 <https://github.com/CWorthy-ocean/ucla-roms/pull/343>`_)

Improvements
~~~~~~~~~~~~

- N/A

Miscellaneous
~~~~~~~~~~~~~

- N/A
