#!/bin/bash
# Note: This script is meant to be run from the root project directory, like this `tools/generate_docs.sh`.
pdoc -f --html -o docs .
