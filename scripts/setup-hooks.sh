#!/bin/bash

# Script to set up git pre-commit hooks for Metheor
# Run this script from the project root directory

set -e

echo "🔧 Setting up pre-commit hooks for Metheor..."

# Check if we're in a git repository
if [ ! -d ".git" ]; then
    echo "❌ Error: This script must be run from the root of a git repository"
    exit 1
fi

# Create the pre-commit hook
cat > .git/hooks/pre-commit << 'EOF'
#!/bin/bash

# Pre-commit hook for Metheor Rust project
# Runs cargo fmt, clippy, and tests before allowing commit

set -e

echo "🔍 Running pre-commit checks..."

# Check if cargo is available
if ! command -v cargo &> /dev/null; then
    echo "❌ cargo is not installed or not in PATH"
    exit 1
fi

echo "📝 Checking code formatting with cargo fmt..."
if ! cargo fmt -- --check; then
    echo "❌ Code formatting issues found!"
    echo "💡 Run 'cargo fmt' to fix formatting issues"
    exit 1
fi
echo "✅ Code formatting is correct"

echo "🔍 Running clippy for linting..."
if ! cargo clippy --all-targets --all-features -- -D warnings; then
    echo "❌ Clippy found issues!"
    echo "💡 Fix the clippy warnings before committing"
    exit 1
fi
echo "✅ Clippy checks passed"

echo "🧪 Running tests..."
if ! cargo test --all-features; then
    echo "❌ Tests failed!"
    echo "💡 Fix failing tests before committing"
    exit 1
fi
echo "✅ All tests passed"

echo "🎉 All pre-commit checks passed! Proceeding with commit..."
EOF

# Make the hook executable
chmod +x .git/hooks/pre-commit

echo "✅ Pre-commit hook installed successfully!"
echo ""
echo "The hook will now run automatically before each commit and will:"
echo "  - Check code formatting with 'cargo fmt'"
echo "  - Run linting with 'cargo clippy'"
echo "  - Run all tests with 'cargo test'"
echo ""
echo "If any check fails, the commit will be blocked until issues are resolved."