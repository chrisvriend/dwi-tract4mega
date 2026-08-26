# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------

project = "TractoPrep"
copyright = "2026, Chris Vriend - Amsterdam UMC"
author = "Chris Vriend - Amsterdam UMC"

import subprocess

try:
    release = subprocess.check_output(
        ["git", "describe", "--tags", "--abbrev=0"], text=True
    ).strip().lstrip("v")
except Exception:
    release = "v1.0.6"

version = release

def replace_release_placeholder(app, docname, source):
    source[0] = source[0].replace("{{RELEASE_TAG}}", release)


def setup(app):
    app.connect("source-read", replace_release_placeholder)

# -- General configuration ----------------------------------------------------

extensions = [
    "myst_parser",
    "sphinx_design",
    "sphinx_substitution_extensions",
    "sphinx.ext.todo",
]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

myst_enable_extensions = [
    "colon_fence",
    "deflist",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Options for HTML output ---------------------------------------------------

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_theme_options = {
    "collapse_navigation": False,
    "navigation_depth": 3,
}

todo_include_todos = True
