.. _{{ objname }}:

{{ fullname | escape | underline }}

.. currentmodule:: {{ module }}
.. autoclass:: {{ objname }}

{% block attributes %}
{% if attributes %}

Attributes table
~~~~~~~~~~~~~~~~
.. autosummary::

{% for item in attributes %}
    ~{{ fullname }}.{{ item }}

{% endfor %}
{% endif %}
{% endblock %}

{% block methods %}
{% if methods %}

Methods table
~~~~~~~~~~~~~
.. autosummary::

{% for item in methods %}
    {%- if item not in inherited_members %}

    ~{{ fullname }}.{{ item }}
    {% endif %}

{% endfor %}
{% endif %}
{% endblock %}

{% block attributes_documentation %}
{% if attributes %}

Attributes
~~~~~~~~~~
{% for item in attributes %}
{{ item | escape | underline(line='^') }}
.. autoattribute:: {{ [objname, item] | join(".") }}
{% endfor %}
{% endif %}
{% endblock %}

{% block methods_documentation %}
{% if methods %}

Methods
~~~~~~~

{% for item in methods %}
{%- if item not in inherited_members %}
.. automethod:: {{ [objname, item] | join(".") }}
{% endif %}
{% endfor %}

Inherited Methods
+++++++++++++++++
{% for item in methods %}
{% if item in inherited_members %}
.. automethod:: {{ [objname, item] | join(".") }}
{% endif %}
{% endfor %}

{% endif %}
{% endblock %}