#!/usr/bin/env bash
# Publish codes and documentation
# Usage: ./publish.sh [local|github] [commit message]

set -euo pipefail

# Default mode is 'local'
MODE=${1:-local}
COMMIT_MSG=${2:-"Update docs"}

echo "🚀 Building CRO documentation in mode: $MODE"

# Build docs
./docsrc_build.sh "$MODE"

if [ "$MODE" = "local" ]; then
    echo "🌐 Serving documentation locally..."
    make serve
else
    echo "📦 Committing and pushing documentation to GitHub..."
    git add docs docsrc pyCRO mCRO
    git commit -m "$COMMIT_MSG"
    git push
fi

echo "✅ Done!"
