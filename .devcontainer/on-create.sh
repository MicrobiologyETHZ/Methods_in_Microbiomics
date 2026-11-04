#!/usr/bin/env bash
# Runs once when the Codespace is created, from the repository root.
set -eu

micromamba env create -y -f .devcontainer/environment.yml

# Every terminal: activate the tools, and start in the workshop folder. (VS Code's
# terminal.integrated.cwd setting doesn't apply to the first terminal Codespaces opens.)
cat >> ~/.bashrc <<BASHRC

micromamba activate gene-catalogs
if [ "\$PWD" = "$PWD" ]; then cd workshop/gene-catalogs; fi
BASHRC
