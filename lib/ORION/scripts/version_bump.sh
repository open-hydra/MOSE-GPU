#!/bin/bash

# Set the file containing the version
VERSION_FILE="scripts/version.txt"

# Function to increment version numbers
increment_version() {
  local version=$1
  local part=$2

  IFS='.' read -r -a parts <<< "$version"

  if [[ "$part" == "major" ]]; then
    parts[0]=$((parts[0] + 1))
    parts[1]=0
    parts[2]=0
  elif [[ "$part" == "minor" ]]; then
    parts[1]=$((parts[1] + 1))
    parts[2]=0
  elif [[ "$part" == "patch" ]]; then
    parts[2]=$((parts[2] + 1))
  else
    echo "Invalid part: $part. Must be 'major', 'minor', or 'patch'."
    exit 1
  fi

  echo "${parts[0]}.${parts[1]}.${parts[2]}"
}

# Ensure the version file exists
if [ ! -f "$VERSION_FILE" ]; then
  echo "Version file $VERSION_FILE not found!"
  exit 1
fi

# Read the current version
current_version=$(cat $VERSION_FILE)
echo "Current version: $current_version"

# Determine which part to increment
if [[ "$1" == "--major" ]]; then
  new_version=$(increment_version "$current_version" "major")
elif [[ "$1" == "--minor" ]]; then
  new_version=$(increment_version "$current_version" "minor")
elif [[ "$1" == "--patch" ]]; then
  new_version=$(increment_version "$current_version" "patch")
else
  echo "Usage: $0 --major|--minor|--patch"
  exit 1
fi

echo "New version: $new_version"

# Update the version file
echo $new_version > $VERSION_FILE

# Commit the version change
git add $VERSION_FILE
git commit -m "Bump version to $new_version"

# Delete all local tags
echo "Deleting all local tags..."
git tag | xargs -r git tag -d

# Delete all remote tags
echo "Deleting all remote tags..."
git ls-remote --tags origin | awk '/refs\/tags\// {print $2}' | while read -r tag; do
  tag_name=${tag#refs/tags/}
  git push origin ":refs/tags/$tag_name"
done

# Create a new Git tag
git tag "v$new_version"

# Push the changes and the new tag
git push origin HEAD:main
git push origin "v$new_version"

echo "Version bumped to $new_version and pushed with tag v$new_version"
