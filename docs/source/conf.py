# Configuration file for the Sphinx documentation builder.

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../../src")))

# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'AlleleTools'
copyright = '2026, Nicolás Mendoza Mejía'
author = 'Nicolás Mendoza Mejía'
release = '0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

# Optional runtime dependencies:
autodoc_mock_imports = ["swifter"]

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

html_favicon = '_static/logo.svg'

# GitHub integration - makes "View page source" link to GitHub
html_context = {
    "display_github": True,
    "github_user": "nmendozam",
    "github_repo": "alleleTools",
    "github_version": "main",
    "conf_py_path": "/docs/source/",
}
