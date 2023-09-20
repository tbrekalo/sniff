#!/usr/bin/env bash
set -eou pipefail

# taken from: https://github.com/foonathan/docker/blob/main/common/install-cmake.sh
#
MIRROR_URL="https://github.com/Kitware/CMake/releases/download/v3.21.4/"
DOWNLOAD_X86="cmake-3.21.4-linux-x86_64.sh"
DOWNLOAD_ARM="cmake-3.21.4-linux-aarch64.sh"
DOWNLOAD_FILE="cmake.sh"

if [ $(uname -m) = "x86_64" ]; then
    DOWNLOAD=$DOWNLOAD_X86
elif [ $(uname -m) = "aarch64" ]; then
    DOWNLOAD=$DOWNLOAD_ARM
else
    echo "unknown architecture" >/dev/stderr
    exit 1
fi

# Download
echo "Downloading ${DOWNLOAD}"
wget -nv -O "${DOWNLOAD_FILE}" "${MIRROR_URL}/${DOWNLOAD}"

# Install
echo "Installing CMakee"
bash "${DOWNLOAD_FILE}" --skip-license --prefix=/usr/local --exclude-subdir

# Cleanup
rm "${DOWNLOAD_FILE}"
