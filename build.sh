#!/usr/bin/env bash
# macOS visualizer: install CMake, GLFW and GLM (e.g. brew install cmake glfw glm).
# Headless core/tests on other platforms: ./build.sh -DBUILD_VISUALIZER=OFF
set -euo pipefail

project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
build_dir="${ORBITAL_BUILD_DIR:-${project_dir}/build}"
if ! command -v cmake >/dev/null 2>&1; then
    echo "CMake is required (on macOS: brew install cmake)." >&2
    exit 1
fi

cmake -S "$project_dir" -B "$build_dir" -DCMAKE_BUILD_TYPE=Release "$@"
cmake --build "$build_dir" --parallel
(cd "$build_dir" && ctest --output-on-failure)
if [[ -x "$build_dir/quantum_orbitals" ]]; then
    echo "Run: $build_dir/quantum_orbitals"
    echo "Controls and scientific conventions: $build_dir/quantum_orbitals --help"
fi
