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
import os
import sys

sys.path.insert(0, os.path.abspath('..'))

# -- Project information -----------------------------------------------------

project = 'CLIVAR Working Group Community Recharge Oscillator (CRO) project'
author  = 'CLIVAR CRO Team (lead by Soong-Ki Kim and Sen Zhao)'
copyright = '2025'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = [
    '_build',
    'Thumbs.db',
    '.DS_Store',
    'notebooks/test-*.ipynb',  # notebooks with filenames starting with "test" won't be rendered
]

extensions = [
    'nbsphinx',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx_copybutton',
    'sphinx_design',
    'sphinxcontrib.bibtex',
    # 'sphinx.ext.mathjax',
    # 'sphinxcontrib.bibtex',
    # 'sphinxcontrib.rsvgconverter',
    # 'sphinx_copybutton',
    # 'sphinx_gallery.load_style',
]
nbsphinx_allow_errors = True
nbsphinx_codecell_lexer = 'python3'

html_logo = 'CRO_logo.png'
html_favicon = 'CRO_logo.png'


# Point to your .bib file(s)
bibtex_bibfiles = ["refs.bib"]

# Optional: APA style
bibtex_default_style = "apa"   # APA via citeproc-py
bibtex_reference_style = "author_year"

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = 'sphinx_rtd_theme'
# html_theme = 'press'
# html_theme = 'pydata_sphinx_theme'
html_theme = 'sphinx_book_theme'
# html_theme = 'sphinx_documatt_theme'
# html_theme = 'furo'
# html_theme = 'sphinx_material'
# html_theme = 'bootstrap'
# html_theme = "sphinxawesome_theme"
html_theme_options = {
    'repository_url': 'https://github.com/senclimate/CRO',
    'use_edit_page_button': True,
    'use_repository_button': True,
    'use_issues_button': True,
    'use_fullscreen_button': True,
    'extra_footer': '<em>CLIVAR Working Group Community Recharge Oscillator (CRO) project</em>',
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".

html_static_path = ['_static']
html_css_files = ['style.css']

def setup(app):
    app.add_css_file('theme_overrides.css')
