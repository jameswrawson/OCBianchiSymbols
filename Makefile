FILES = $(wildcard *.sage)
OBJECTS = $(FILES:.sage=.py)

all: $(OBJECTS)

%.py: %.sage
	sage --preparse $<
	mv $(addsuffix .py,$<) $@
