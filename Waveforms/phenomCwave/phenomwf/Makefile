testphenomwf_DEPS = 	phenomwf.c \
			phenomwf.h

default: testphenomwf

clean:
	rm -f testphenomwf 

testphenomwf: testphenomwf.c ${testphenomwf_DEPS}
	gcc -o testphenomwf -I. -lm testphenomwf.c phenomwf.c
