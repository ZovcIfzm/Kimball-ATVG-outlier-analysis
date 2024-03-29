# Signifies our desired python version
# Makefile macros (or variables) are defined a little bit differently than traditional bash, keep in mind that in the Makefile there's top-level Makefile-only syntax, and everything else is bash script syntax.
PYTHON = python3

# .PHONY defines parts of the makefile that are not dependant on any specific file
# This is most often used to store functions
.PHONY = help

# Defines the default target that `make` will to try to make, or in the case of a phony target, execute the specified commands
# This target is executed whenever we just type `make`
.DEFAULT_GOAL = help

# The @ makes sure that the command itself isn't echoed in the terminal
help:
	@echo "---------------HELP-----------------"
	@echo "To run the code type make run"
	@echo "To run without plots type make no_plot"
	@echo "------------------------------------"

# The ${} notation is specific to the make syntax and is very similar to bash's $() 
# This function uses pytest to test our source files
run:
	pip install -r requirements.txt
	${PYTHON} clustering_outliers.py --counts_file results/count_output.txt --show_plots 1

no_plot:
	pip install -r requirements.txt
	${PYTHON} clustering_outliers.py --counts_file results/count_output.txt --show_plots 0
