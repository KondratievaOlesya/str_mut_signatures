# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'str_mut_signatures'
copyright = '2026, Olesia Kondrateva'
author = 'Olesia Kondrateva'
release = '2.0.0'
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",      # NumPy/Google docstrings
    "sphinx.ext.autosummary",   # generates API stub pages
    "sphinx.ext.viewcode",      # link to source
    "myst_parser",              # allow Markdown pages
    "sphinxarg.ext",
]

autosummary_generate = True

napoleon_numpy_docstring = True
napoleon_google_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_param = True
napoleon_use_rtype = True

autodoc_member_order = "bysource"
autodoc_typehints = "signature"  # show type hints in signature, keep docstrings clean

html_theme = "sphinx_rtd_theme"

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output
html_static_path = ['_static']
