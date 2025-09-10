# Root Makefile for CRO docs

SPHINXBUILD ?= sphinx-build
SPHINXOPTS  ?=

# Source directories
PY_SRC = docsrc/python
MAT_SRC = docsrc/matlab

# Output directories
PY_OUT = docs/python
MAT_OUT = docs/matlab

.PHONY: help python matlab docs serve clean

help:
	@echo "make python        - build Python docs"
	@echo "make matlab        - build Matlab docs"
	@echo "make docs          - build both"
	@echo "make serve         - serve locally at http://localhost:8000"
	@echo "make clean         - remove all built docs"

python:
	@echo "🟢 Building Python docs..."
	$(SPHINXBUILD) -b html $(PY_SRC) $(PY_OUT) $(SPHINXOPTS)

matlab:
	@echo "🟢 Building Matlab docs..."
	$(SPHINXBUILD) -b html $(MAT_SRC) $(MAT_OUT) $(SPHINXOPTS)

docs: python matlab

serve:
	@echo "🌐 Serving docs at http://localhost:8000"
	cd docs && python -m http.server 8000

clean:
	@echo "🧹 Cleaning build directories..."
	rm -rf $(PY_OUT)/* $(MAT_OUT)/*
	rm -f docs/index.html docs/switcher.json
