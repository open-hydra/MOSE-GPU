# ORION

## Versioning

This project follows [Semantic Versioning](https://semver.org/) (SemVer) for managing version numbers. Semantic Versioning is a widely adopted versioning scheme that conveys meaning about the underlying code and what has been modified.

### Version Format

The version number is structured as follows:

MAJOR.MINOR.PATCH

- **MAJOR**: Incremented when there are incompatible API changes.
- **MINOR**: Incremented when functionality is added in a backwards-compatible manner.
- **PATCH**: Incremented when backwards-compatible bug fixes are made.

For example, a version number of `1.2.3` indicates:
- `1`: The first major version, with potential breaking changes since version `0.x.x`.
- `2`: The second minor update, adding new features without breaking existing ones.
- `3`: The third patch, fixing bugs in a backwards-compatible manner.

### How Versioning is Managed

Version numbers are automatically updated running the following command AFTER committing the code changes:

```bash
./scripts/version_bump.sh --major|--minor|--patch
```

Pushing is not required. The updated version is tagged in the repository and released automatically.