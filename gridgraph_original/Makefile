FC = gfortran
FFLAGS = -std=legacy

SRC = ggraph1.f
PATCHED = ggraph1_patched.f
BIN = ggraph

.PHONY: all run clean

all: $(BIN)

# Patch the source before compiling:
#   1. Rename rand() calls/definition/return to schrage_rand to avoid
#      conflict with gfortran's built-in rand() intrinsic.
#   2. Replace etime(it) with 0.0 (timing call not portable to gfortran).
$(PATCHED): $(SRC)
	perl -pe \
	  's/([^a-z_])rand\(/$${1}schrage_rand(/g; \
	   s/\brand=/schrage_rand=/; \
	   s/real rand$$/real schrage_rand/; \
	   s/etime\(it\)/0.0/g' \
	  $(SRC) > $(PATCHED)

$(BIN): $(PATCHED)
	$(FC) $(FFLAGS) -o $(BIN) $(PATCHED)

# Usage: make run ARGS="h w MAXCAP MAXCOST SEED"
# Example: make run ARGS="8 8 10000 1000 270001"
run: $(BIN)
	echo "$(ARGS)" | ./$(BIN)

clean:
	rm -f $(BIN) $(PATCHED)
