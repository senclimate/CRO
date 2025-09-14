#!/usr/bin/env bash
# Build and prepare CRO docs for local or GitHub Pages
# Usage: ./docsrc_build.sh [local|publish]

set -e

MODE=${1:-local}  # default to local if not specified

echo "🚀 Building CRO docs in mode: $MODE"

# Clean
make clean

# Build both projects
make docs

# Prepare switcher.json
SWITCHER_FILE=docs/switcher.json
if [ "$MODE" = "local" ]; then
    cat > "$SWITCHER_FILE" <<EOF
[
  {"version": "Python", "url": "/python/"},
  {"version": "Matlab", "url": "/matlab/"}
]
EOF
    echo "🔹 Local switcher.json created (localhost paths)"
else
    cat > "$SWITCHER_FILE" <<EOF
[
  {"version": "Python", "url": "https://senclimate.github.io/CRO/python/"},
  {"version": "Matlab", "url": "https://senclimate.github.io/CRO/matlab/"}
]
EOF
    echo "🔹 Production switcher.json created (GitHub Pages paths)"
fi

# Create root landing page (redirect to Python docs by default)
cat > docs/index.html <<EOF
<!DOCTYPE html>
<html>
<head>
  <meta charset="UTF-8">
  <title>CRO Documentation</title>
  <meta http-equiv="refresh" content="0; url=python/">
</head>
<body>
  <p>Redirecting to <a href="python/">Python Documentation</a>...</p>
</body>
</html>
EOF
echo "🔹 Root index.html created (Python docs default)"

# GitHub Pages preparation
touch docs/.nojekyll
git add docs

echo "✅ Docs built and ready at docs/ (mode=$MODE)"
