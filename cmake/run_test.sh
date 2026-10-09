#!/bin/sh
# Run one CTest command with a private TMPDIR. Many tests write fixed names
# under the temporary directory, so tests running concurrently (ctest -j, or
# make -j for a test goal) would otherwise overwrite each other's files.
# Set TUBULARFLOWIGA_KEEP_TEST_OUTPUT to keep the directory after the test.
TMPDIR=$(mktemp -d "${TMPDIR:-/tmp}/tfi-test.XXXXXX") || exit 1
export TMPDIR
sh -c "$1"
status=$?
if [ -n "$TUBULARFLOWIGA_KEEP_TEST_OUTPUT" ]; then
	echo "test output kept in $TMPDIR"
else
	rm -rf "$TMPDIR"
fi
exit $status
