
# Simple makefile for dev.

all:
	# Use the local build wrapper to automate writing the report log to stdout.
	./build.sh

clean-build:
	rm -fR _build


clean-files:
	@find . -name ".pytest_cache" -exec rm -rf {} \;
	@find . -name "*undo-tree*" -exec rm -f {} \;


clean: clean-build clean-files
