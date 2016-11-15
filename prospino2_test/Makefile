prospino_dir = $(shell git rev-parse --show-toplevel)/vendor/prospino2

src = $(wildcard pros_*.f90)
bin = $(addsuffix .run, $(basename $(src)))

all: $(bin)

include $(prospino_dir)/Makefile
DIRECT = $(prospino_dir)/
FILES  = ${FILES_DIRECT} ${FILES_MACROS} ${FILES_MATRIX} \
         ${FILES_INTEGS} ${FILES_STRONG} ${FILES_INTERF} 
.SECONDARY: ${FILES}

$(FILES):
	FC=gfortran make -C $(prospino_dir) $@

%.run : %.f90 $(FILES)
	gfortran -I$(prospino_dir) $^ -o $@
