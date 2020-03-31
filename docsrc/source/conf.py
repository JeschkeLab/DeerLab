# Configuration file for the Sphinx documentation builder.
# -- Project information -----------------------------------------------------

import sys
import os
import re

sys.path.append(os.path.abspath('.'))

# Import Read-the-Docs (RTD) theme
from sphinx.locale import _
from sphinx_rtd_theme import __version__

# Project details
project = 'DeerLab'
copyright = '2019, Luis Fabregás Ibáñez, Stefan Stoll, and others'
author = 'Fabregas Ibanez'
language = 'en'

# Add sphinx extensions
extensions = [
    'sphinx.ext.intersphinx',
    'sphinxcontrib.matlab',
    'sphinx.ext.viewcode',
    'sphinxcontrib.httpdomain',
    'sphinx.ext.imgmath',
    'button'
]

# Render Latex math equations as svg instead of rendering with JavaScript
imgmath_image_format = 'svg'
imgmath_dvisvgm = 'dvisvgm'


# Setup template stuff
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

# Integrate version control system
# -------------------------------------------------------------
html_context = {
    "display_github": False, # Integrate GitHub
    "github_user": "luisfabib", # Username
    "github_repo": "deerlab", # Repo name
    "github_version": "master", # Version
    "conf_py_path": "/source/", # Path in the checkout to the docs root
}


# Read-the-Docs options configuration
# -------------------------------------------------------------
html_theme_options = {
    'sticky_navigation': False,
    'titles_only': True,
    'logo_only': True,
    'navigation_depth':0,
    #'style_nav_header_background': '#2a7bf8'
}
html_copy_source = False
html_theme_path = ["../.."]
html_logo = "demo/static/logo-wordmark-light.svg"
html_show_sourcelink = True
html_favicon = 'favicon.ico'
html_static_path = ['_static']

# Extensions to theme docs
def setup(app):
    from sphinx.domains.python import PyField
    from sphinx.util.docfields import Field
    app.add_css_file('/source/_static/custom.css')
    app.add_stylesheet('/source/_static/theme_override.css')
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

# These folders are copied to the documentation's HTML output
html_static_path = ['_static']

# Add path to custom CSS file to overwrite some of the default CSS settings
html_css_files = [
    'custom.css',
    'theme_override.css'
]	

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# Default role
default_role = 'math'  # with this, :math:`\psi` can be written simply as `\psi`


# -- Options for HTML output -------------------------------------------------
html_theme = "sphinx_rtd_theme"
html_theme_path = ["_themes", ]
html_static_path = ['_static']
highlight_language = 'matlab'
primary_domain = 'mat'
html_logo = '../source/logo.png'


# Design pygments patch for MATLAB true code highlighting
# --------------------------------------------------------

# Import the pygments python library
from pygments.style import Style
from pygments.token import Text, Other, Comment, Whitespace
from pygments.style import Style
from pygments.token import Keyword, Name, Comment, String, Error, Number, Operator, Generic

# Define custom style for the MATLAB highlighting
class MyFancyStyle(Style):
    background_color = "#f7f7f7"
    default_style = ""
    styles = {
        Comment:                'italic #22924a',
        Keyword:                '#0000FF',
        Name:                   '#000',
        Name.Function:          '#000',
        Name.Class:             '#000',
        String:                 '#a022f0'
    }
	
# Create patch for applying the style 	
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

