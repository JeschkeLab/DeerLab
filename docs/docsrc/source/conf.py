# Configuration file for the Sphinx documentation builder.
#import os
#import sys
#sys.path.insert(0, os.path.abspath('./venv/Lib/site-packages/sphinxcontrib_matlabdomain-0.8.0.dist-info'))

# -- Project information -----------------------------------------------------

import sys
import os
import re

sys.path.append(os.path.abspath('..'))

#Import RTD theme
from sphinx.locale import _
from sphinx_rtd_theme import __version__

#Project details
project = 'DeerAnalysis2'
copyright = '2019, Luis Fabregas Ibanez'
author = 'Fabregas Ibanez'
language = 'en'
release = '2019'

#Add sphinx extensions
extensions = [
    'sphinx.ext.intersphinx',
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinxcontrib.matlab',
    'sphinxcontrib.httpdomain',
    'sphinx.ext.imgmath',
]

#Render Latex math equations as svg instead of rendering with JavaScript
imgmath_image_format = 'svg'
imgmath_dvisvgm = 'dvisvgm'


#Setup template stuff
templates_path = ['_templates']
source_suffix = '.rst'
exclude_patterns = []
master_doc = 'index'
suppress_warnings = ['image.nonlocal_uri']
pygments_style = 'default'
intersphinx_mapping = {
    'rtd': ('https://docs.readthedocs.io/en/latest/', None),
    'sphinx': ('http://www.sphinx-doc.org/en/stable/', None),
}
html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'logo_only': True
}
html_theme_path = ["../.."]
html_logo = "demo/static/logo-wordmark-light.svg"
html_show_sourcelink = True
htmlhelp_basename = slug
latex_documents = [
  ('index', '{0}.tex'.format(slug), project, author, 'manual'),
]
man_pages = [
    ('index', slug, project, [author], 1)
]
texinfo_documents = [
  ('index', slug, project, author, slug, project, 'Miscellaneous'),
]

# Extensions to theme docs
def setup(app):
    from sphinx.domains.python import PyField
    from sphinx.util.docfields import Field

    app.add_object_type(
        'confval',
        'confval',
        objname='configuration value',
        indextemplate='pair: %s; configuration value',
        doc_field_types=[
            PyField(
                'type',
                label=_('Type'),
                has_arg=False,
                names=('type',),
                bodyrolename='class'
            ),
            Field(
                'default',
                label=_('Default'),
                has_arg=False,
                names=('default',),
            ),
        ]
    )

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
html_theme = "sphinx_rtd_theme"
html_theme_path = ["_themes", ]
html_static_path = ['_static']
highlight_language = 'matlab'
primary_domain = 'mat'
html_logo = '../source/logo.png'


#Design pygments patch for MATLAB true code highlighting
from pygments.style import Style
from pygments.token import Text, Other, Comment, Whitespace
from pygments.style import Style
from pygments.token import Keyword, Name, Comment, String, Error, Number, Operator, Generic

class MyFancyStyle(Style):
    background_color = "#ECEBEB"
    default_style = ""
    styles = {
        Comment:                'italic #4B7F4B',
        Keyword:                '#000',
        Name:                   '#000',
        Name.Function:          '#000',
        Name.Class:             '#000',
        String:                 '#CC00FF'
    }
def pygments_monkeypatch_style(mod_name, cls):
    import sys
    import pygments.styles
    cls_name = cls.__name__
    mod = type(__import__("os"))(mod_name)
    setattr(mod, cls_name, cls)
    setattr(pygments.styles, mod_name, mod)
    sys.modules["pygments.styles." + mod_name] = mod
    from pygments.styles import STYLE_MAP
    STYLE_MAP[mod_name] = mod_name + "::" + cls_name
pygments_monkeypatch_style("my_fancy_style", MyFancyStyle)
pygments_style = "my_fancy_style"

