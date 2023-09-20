#!/bin/bash
set -eou pipefail

apt-get update && apt-get install \
  -y --no-install-recommends git ninja-build pkg-config python3-pip
pip3 install conan==2.0

conan profile detect -f
