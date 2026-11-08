#!/bin/bash
# Fetch AutoQ-Para at the pinned commit, apply our build fixes, add our models (autoq_models.cpp) and build.
# Needs git, cmake >= 3.25, g++ and libgmp-dev. Catch2 is not needed (its tests are not built).
# Usage: bash experiments/autoq_para/fetch_and_build.sh [DIR]   (default third_party_build/AutoQ-Para)
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
DIR=${1:-"$HERE/../../third_party_build/AutoQ-Para"}
COMMIT=e3a4718
[ -d "$DIR/.git" ] || git clone https://github.com/VeriFIT/AutoQ-Para "$DIR"
cd "$DIR"
git checkout -q "$COMMIT"
git apply --check "$HERE/build_fixes.patch" 2>/dev/null && git apply "$HERE/build_fixes.patch"
cp "$HERE/autoq_models.cpp" src/
cmake -S src -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target autoq_models -j "$(nproc)"
echo "built $DIR/build/autoq_models"
