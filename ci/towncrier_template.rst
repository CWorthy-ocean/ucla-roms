{#-
  towncrier template for docs/releases/<tag>.rst (see towncrier.toml).
  It writes the ".. _<tag>:" anchor and the title, then a heading for each
  fragment type that has fragments and one "- " bullet per fragment, laid out
  as in the releases already in docs/releases/. The PR link goes at the end of
  the fragment's first paragraph. Whatever follows it (sub-bullets, further
  paragraphs) comes after a blank line, which RST needs before a nested list,
  and a bullet that had such a tail is followed by a blank line too.
  A fragment written as a bullet ("- Fixed X") loses its marker, and the two
  columns of indentation that went with it, so it does not render as a list
  nested in the bullet. A first paragraph ending in "::" opens a literal block,
  so its link goes in a paragraph of its own after the block. (towncrier pads
  such a fragment with a line holding one space; rstrip drops it.)
-#}
.. _{{ versiondata.version }}:

{{ versiondata.version }}
{{ top_underline * versiondata.version | length }}
{% for section, _ in sections.items() %}
{% if sections[section] %}
{% for category, val in definitions.items() if category in sections[section] %}

{{ definitions[category]["name"] }}
{{ underlines[0] * definitions[category]["name"] | length }}

{% set state = namespace(after_tail=false) %}
{% for text, values in sections[section][category].items() %}
{% set entry = namespace(lead=[], tail=[], blank=false, bullet=false) %}
{% for line in text.rstrip().split("\n") %}
{% if loop.first %}
{% set entry.bullet = line[:2] in ("- ", "* ") %}
{% set entry.lead = [line[2:] if entry.bullet else line] %}
{% else %}
{% set line = line[2:] if entry.bullet and line[:4] == "    " else line %}
{% if entry.tail or (entry.blank and line.strip()) or line.lstrip()[:2] in ("- ", "* ") %}
{% set entry.tail = entry.tail + [line] %}
{% elif not line.strip() %}
{% set entry.blank = true %}
{% else %}
{% set entry.lead = entry.lead + [line] %}
{% endif %}
{% endif %}
{% endfor %}
{% set link = "(" + values | join(", ") + ")" if values else "" %}
{% set literal = entry.lead[-1].rstrip()[-2:] == "::" %}
{% if state.after_tail %}

{% endif %}
- {{ entry.lead | join("\n") }}{% if link and not literal %} {{ link }}{% endif %}

{% if entry.tail %}

{{ entry.tail | join("\n") }}
{% endif %}
{% if link and literal %}

  {{ link }}
{% endif %}
{% set state.after_tail = (entry.tail | length > 0) or (link and literal) %}
{% endfor %}
{% endfor %}
{% else %}

No significant changes.
{% endif %}
{% endfor %}
