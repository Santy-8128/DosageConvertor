#!/bin/bash


rm -rf cget/ release-build/ install.log
if [[ ! $(which cget) ]]; then
  (>&2 echo "Error: cget not installed. Please run 'pip install --user cget'.")
  exit -1
fi

echo -e "Installing Dependencies - Libstatgen ..."
cget install -f requirements.txt
rc=$?
if [[ $rc != 0 ]]; then
  (>&2 echo "Failed Installing Dependencies - Libstatgen")
  exit $rc
fi

mkdir release-build
cd release-build/
echo -e "Configuring ..."
cmake -DCMAKE_TOOLCHAIN_FILE=../cget/cget/cget.cmake -DCMAKE_BUILD_TYPE=Release .. || exit -1

echo -e "Building ..."
make || exit -1

echo "Binary created at /release-build/DosageConvertor"

