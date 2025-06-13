#Look for an executable called cffconvert
find_program(cffconvert
             NAMES cffconvert
             DOC "Path to sphinx-build executable")

include(FindPackageHandleStandardArgs)

#Handle standard arguments to find_package like REQUIRED and QUIET
find_package_handle_standard_args(cffconvert
                                  "Failed to find cffconvert executable"
                                  cffconvert)