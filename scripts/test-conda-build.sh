#!/bin/bash
set -e

echo "🧪 Testing Metheor Conda Build System"
echo "======================================"

# Check if rattler-build is available
if ! command -v rattler-build &> /dev/null; then
    echo "❌ rattler-build not found. Installing with pixi..."
    
    # Check if pixi is available
    if ! command -v pixi &> /dev/null; then
        echo "❌ pixi not found. Please install pixi first:"
        echo "   curl -fsSL https://pixi.sh/install.sh | bash"
        exit 1
    fi
    
    pixi global install rattler-build
fi

echo "✅ rattler-build found: $(rattler-build --version)"

# Check if recipe exists
if [ ! -f "recipe.yaml" ]; then
    echo "❌ recipe.yaml not found in current directory"
    exit 1
fi

echo "✅ recipe.yaml found"

# Validate recipe syntax
echo "🔍 Validating recipe syntax..."
if rattler-build build --recipe recipe.yaml --dry-run; then
    echo "✅ Recipe syntax is valid"
else
    echo "❌ Recipe syntax validation failed"
    exit 1
fi

# Test local build for current platform
echo "🔨 Testing local build..."
platform=$(uname -m)
case $(uname -s) in
    Linux*)     target_platform="linux-64";;
    Darwin*)    
        if [ "$platform" = "arm64" ]; then
            target_platform="osx-arm64"
        else
            target_platform="osx-64"
        fi;;
    MINGW*)     target_platform="win-64";;
    *)          target_platform="unknown";;
esac

echo "🎯 Building for platform: $target_platform"

# Create output directory
mkdir -p test-build-output

# Build the package
if rattler-build build \
    --recipe recipe.yaml \
    --target-platform "$target_platform" \
    --channel conda-forge \
    --output-dir ./test-build-output; then
    
    echo "✅ Local build successful!"
    
    # List built packages
    echo "📦 Built packages:"
    find test-build-output -name "*.conda" -o -name "*.tar.bz2" | while read package; do
        echo "   - $(basename "$package")"
    done
    
else
    echo "❌ Local build failed"
    exit 1
fi

echo ""
echo "🎉 All tests passed! The conda build system is ready."
echo ""
echo "Next steps:"
echo "1. Set up ANACONDA_API_TOKEN secret in GitHub (see CONDA_BUILD_SETUP.md)"
echo "2. Create a new git tag: git tag v0.1.9 && git push origin v0.1.9"
echo "3. Monitor the build at: https://github.com/dohlee/metheor/actions"

# Cleanup
echo "🧹 Cleaning up test build artifacts..."
rm -rf test-build-output