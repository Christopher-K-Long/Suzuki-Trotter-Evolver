CITATION_VERSION=$(cat "../CITATION.cff" | perl -nle'print $& while m{(?<=^version:\s).*}g')

CHANGELOG_VERSION=$(cat "../ChangeLog.md" | perl -nle'print $& while m{(?<=^##\sRelease\s).*}g' | head -1)

README_VERSION_1=$(cat "../README.md" | perl -nle'print $& while m{(?<=\[`).*(?=`\]\(ChangeLog\.md#release-)}g')
README_VERSION_2=$(cat "../README.md" | perl -nle'print $& while m{(?<=`\]\(ChangeLog\.md#release-)[0-9]*(?=\))}g')

CMAKE_PROJECT_VERSION=$(cat "../CMakeLists.txt" | perl -nle'print $& while m{(?<=project\(Suzuki-Trotter-Evolver VERSION ).*}g')

VCPKG_VERSION=$(cat "../vcpkg.json" | perl -nle'print $& while m{(?<="version":\s").*(?=")}g')

PASSING=true

if [[ "$CMAKE_PROJECT_VERSION" != "$CITATION_VERSION" ]]
then
echo "version in CITATION.cff (${CITATION_VERSION}) does not match cmake project version (${CMAKE_PROJECT_VERSION})"
PASSING=false
fi

if [[ "$CMAKE_PROJECT_VERSION" != "$CHANGELOG_VERSION" ]]
then
echo "version in ChangeLog.md (${CHANGELOG_VERSION}) does not match cmake project version (${CMAKE_PROJECT_VERSION})"
PASSING=false
fi

if [[ "$CMAKE_PROJECT_VERSION" != "$README_VERSION_1" ]]
then
echo "version in README.md (${README_VERSION_1}) does not match cmake project version (${CMAKE_PROJECT_VERSION})"
PASSING=false
fi

if [[ "$CMAKE_PROJECT_VERSION" != "$VCPKG_VERSION" ]]
then
echo "version in vcpkg.json (${VCPKG_VERSION}) does not match cmake project version (${CMAKE_PROJECT_VERSION})"
PASSING=false
fi

FORMATTED_CMAKE_PROJECT_VERSION=$(echo "$CMAKE_PROJECT_VERSION" | sed 's/[^0-9]*//g')

if [[ "$FORMATTED_CMAKE_PROJECT_VERSION" != "$README_VERSION_2" ]]
then
echo "version in README.md (${README_VERSION_2}) does not match cmake project version (${CMAKE_PROJECT_VERSION})"
PASSING=false
fi

if [[ "$PASSING" == false ]]
then
    exit 1
fi

echo "Passed version consistency test: all versions match."