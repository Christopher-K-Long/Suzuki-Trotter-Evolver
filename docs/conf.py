# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Suzuki-Trotter-Evolver'
copyright = '2025, Christopher K. Long'
author = 'Christopher K. Long'

# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["myst_parser",
              "breathe",
              "exhale",
              "sphinx.ext.mathjax",
              "sphinx_tabs.tabs"]
myst_enable_extensions = ["dollarmath"]

myst_heading_anchors = 5
suppress_warnings = ["myst.header"]
breathe_default_project = "Suzuki-Trotter-Evolver"
exhale_args = {
    # These arguments are required
    "containmentFolder":     "reference",
    "rootFileName":          "index.rst",
    "doxygenStripFromPath":  "..",
    # Heavily encouraged optional argument (see docs)
    "rootFileTitle":         "API Reference",
    # Suggested optional arguments
    "contentsDirectives" : False,
    "createTreeView":        False,
    # TIP: if using the sphinx-bootstrap-theme, you need
    # "treeViewIsBootstrap": True,
    # "exhaleExecutesDoxygen": True,
    # "exhaleDoxygenStdin":    "INPUT = ../include"
}
primary_domain = 'cpp'
highlight_language = 'cpp'
mathjax3_config = {'loader': {'load': ['[tex]/mathtools', '[tex]/physics']},
                   'tex': {'packages': {'[+]': ['mathtools', 'physics']}},}

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'furo'
html_static_path = ['_static']


import subprocess, os, cffconvert

DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
cffconvert_command = os.path.join(cffconvert.__file__.split("lib")[0], "bin", "cffconvert")
subprocess.check_call(f"cd {DIR}; {cffconvert_command} -f apalike > docs/citation/citation_files/citation.txt", shell=True)
subprocess.check_call(f"cd {DIR}; {cffconvert_command} -f bibtex > docs/citation/citation_files/citation.bib", shell=True)
subprocess.check_call(f"cd {DIR}; {cffconvert_command} -f ris > docs/citation/citation_files/citation.ris", shell=True)
subprocess.check_call(f"cd {DIR}; {cffconvert_command} -f codemeta > docs/citation/citation_files/citation_codemeta.json", shell=True)
subprocess.check_call(f"cd {DIR}; {cffconvert_command} -f endnote > docs/citation/citation_files/citation.enw", shell=True)

def configureDoxyfile(input_dir, output_dir):
    with open('Doxyfile.in', 'r') as file :
        filedata = file.read()

    filedata = filedata.replace('@DOXYGEN_INPUT_DIR@', input_dir)
    filedata = filedata.replace('@DOXYGEN_OUTPUT_DIR@', output_dir)

    with open('Doxyfile', 'w') as file:
        file.write(filedata)

# Check if we're running on Read the Docs' servers
read_the_docs_build = os.environ.get('READTHEDOCS', None) == 'True'

if read_the_docs_build:
    breathe_projects = {}
    input_dir = '../include/Suzuki-Trotter-Evolver'
    output_dir = 'build'
    configureDoxyfile(input_dir, output_dir)
    subprocess.call('doxygen', shell=True)
    breathe_projects['Suzuki-Trotter-Evolver'] = output_dir + '/xml'