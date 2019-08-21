# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
#import os
#import sys
#sys.path.insert(0, os.path.abspath('./venv/Lib/site-packages/sphinxcontrib_matlabdomain-0.8.0.dist-info'))

# -- Project information -----------------------------------------------------

import sys
import os
import re

sys.path.append(os.path.abspath('..'))
sys.path.append(os.path.abspath('./demo/'))

from sphinx.locale import _

from sphinx_rtd_theme import __version__


project = u'Read the Docs Sphinx Theme'
slug = re.sub(r'\W+', '-', project.lower())
version = __version__
release = __version__
author = u'Dave Snider, Read the Docs, Inc. & contributors'
copyright = author
language = 'en'

extensions = [
    'sphinx.ext.intersphinx',
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinxcontrib.matlab',
    'sphinxcontrib.httpdomain',
    'sphinx.ext.imgmath',
]
#extensions = [
#    'sphinx.ext.intersphinx',
#    'sphinx.ext.autodoc',
#   # 'sphinx.ext.mathjax',
#    'sphinx.ext.viewcode',
#    'sphinxcontrib.httpdomain',
#   # 'sphinxcontrib.katex',
#    'sphinx.ext.imgmath'
#]



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
project = 'DeerAnalysis2'
copyright = '2019, Luis Fabregas Ibanez'
author = 'Fabregas Ibanez'

# The full version, including alpha/beta/rc tags
release = '2019'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.

import sphinx_rtd_theme


# Add any paths that contain templates here, relative to this directory.
templates_path = [
    '_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"
html_theme_path = ["_themes", ]
# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
highlight_language = 'matlab'
primary_domain = 'mat'
html_logo = '../source/logo.png'



# BEGIN MONKEY-PATCH
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
# END MONKEY-PATCH

imgmath_image_format = 'svg'
imgmath_dvisvgm = 'dvisvgm'
