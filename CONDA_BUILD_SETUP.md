# Automated Conda Build Setup

This guide explains how to set up the automated multi-platform conda build system for Metheor.

## Overview

The automated build system:
- ✅ **Builds for 4 platforms**: `linux-64`, `osx-64`, `osx-arm64`, `win-64`
- ✅ **Uses rattler-build**: Modern, fast alternative to conda-build
- ✅ **Triggers on git tags**: Automatic builds when you create version tags
- ✅ **Uploads to anaconda.org**: Automatically publishes packages
- ✅ **Creates GitHub releases**: Attaches conda packages to releases

## Setup Instructions

### 1. Configure Anaconda API Token

You need to set up an API token for uploading to anaconda.org:

#### Get your Anaconda API token:
```bash
# Login to anaconda.org
anaconda login

# Generate a new API token
anaconda auth --create --name "GitHub-Actions-Metheor" --scopes "api:write api:read"
```

#### Add the token to GitHub Secrets:
1. Go to your GitHub repository: `https://github.com/dohlee/metheor`
2. Navigate to: **Settings** → **Secrets and variables** → **Actions**
3. Click **"New repository secret"**
4. Name: `ANACONDA_API_TOKEN`
5. Value: Paste your token from the previous step
6. Click **"Add secret"**

### 2. Test the Build System

#### Option A: Test with a new tag (recommended)
```bash
# Create and push a new tag
git tag v0.1.9
git push origin v0.1.9
```

#### Option B: Manual trigger
1. Go to **Actions** tab in GitHub
2. Select **"Build Conda Packages"** workflow
3. Click **"Run workflow"** → **"Run workflow"**

### 3. Monitor the Build

After triggering a build:

1. **Check GitHub Actions**: Go to the Actions tab to monitor progress
2. **View build logs**: Click on the running workflow to see detailed logs
3. **Check artifacts**: Completed builds will have downloadable conda packages
4. **Verify upload**: Check https://anaconda.org/dohlee/metheor for new packages

## File Structure

The new build system adds these files:

```
metheor/
├── recipe.yaml                           # Rattler-build recipe (replaces meta.yaml)
├── .github/workflows/conda-build.yml     # Multi-platform build workflow  
├── conda-recipes/                        # Old conda-build files (can be removed)
│   ├── meta.yaml                         # Legacy recipe
│   └── build.sh                          # Legacy build script
└── CONDA_BUILD_SETUP.md                  # This setup guide
```

## Supported Platforms

| Platform | OS | Architecture | Status |
|----------|-------|-------------|--------|
| `linux-64` | Ubuntu 20.04 | x86_64 | ✅ Supported |
| `osx-64` | macOS 13 | Intel | ✅ Supported |
| `osx-arm64` | macOS 14 | Apple Silicon | ✅ Supported |
| `win-64` | Windows | x86_64 | ✅ Supported |

## Build Triggers

The build system triggers on:

- **Git tags** starting with `v*` (e.g., `v0.1.9`) → **Builds + Uploads**
- **Pull requests** to `master`/`dev` branches → **Test builds only**
- **Manual dispatch** → **Test builds only**

## Troubleshooting

### Build fails on specific platform
- Check the Actions logs for the failing platform
- Common issues: missing system dependencies, compilation errors
- The build matrix allows other platforms to continue even if one fails

### Upload fails
- Verify `ANACONDA_API_TOKEN` secret is set correctly
- Check that the token has `api:write` permissions
- Duplicate uploads are automatically skipped

### Version mismatch
- Ensure the git tag matches the version in `recipe.yaml`
- Update the `context.version` field in `recipe.yaml` if needed

## Migration from Old System

The old conda-build system in `conda-recipes/` can be safely removed after verifying the new system works:

```bash
# After successful test, remove old files
rm -rf conda-recipes/
```

## Advanced Configuration

### Custom version from git
The recipe supports automatic versioning from git tags:

```yaml
context:
  version: ${{ git.latest_tag() | replace("v", "") }}
```

### Additional build options
See [rattler-build documentation](https://rattler.build) for advanced configuration options.

## Support

If you encounter issues:
1. Check GitHub Actions logs first
2. Verify all secrets are configured correctly  
3. Test locally with: `rattler-build build --recipe recipe.yaml`