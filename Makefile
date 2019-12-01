include config.mk
CFLAGS := -I$(GMP)/include -I$(MPFR)/include $(CFLAGS)
OBJS = dicho.o gl_gen.o interval.o nc_gen.o newton_root.o usp.o utils.o\
	secante_root.o adaptive.o listz.o polyeval.o median.o frac_rec.o
LIBS = libcrq.a
DIST = Makefile dicho.c gl_gen.c interval.c interval.h crq-impl.h crq_structs.h\
	crq.h nc_gen.c newton_root.c usp.c usp.h utils.c utils.h check.c\
	secante_root.c config.mk COPYING.LIB README adaptive.c filter nc_bench.c\
	listz.c polyeval.c median.c frac_rec.c listz.h
VERSION = 0.0.20061023

all: check.static $(LIBS)

libcrq.a: $(OBJS)
	ar cru libcrq.a $(OBJS)
	ranlib libcrq.a

check: check.o libcrq.a
	gcc $(CFLAGS) $^ -L$(GMP)/lib -L$(MPFR)/lib -lgmp -lmpfr -o $@

exemple: exemple.o libcrq.a
	gcc $(CFLAGS) $^ -L$(GMP)/lib -L$(MPFR)/lib -lgmp -lmpfr -o $@

check.static: check.o libcrq.a
	gcc $(CFLAGS) $^ $(MPFR)/lib/libmpfr.a $(GMP)/lib/libgmp.a -o $@

nc_bench: nc_bench.o libcrq.a
	gcc $(CFLAGS) $^ $(MPFR)/lib/libmpfr.a $(GMP)/lib/libgmp.a -o $@

pi: pi.o libcrq.a
	gcc $(CFLAGS) $^ -L$(GMP)/lib -L$(MPFR)/lib -lgmp -lmpfr -o $@

install: $(LIBS)
	$(INSTALL) --mode=644 crq.h $(DESTDIR)/usr/include
	$(INSTALL) libcrq.a $(DESTDIR)/usr/lib

clean:
	rm -f $(OBJS) $(LIBS) check check.o

dist:	crq-$(VERSION).tar.gz

publish: dist
	scp crq-$(VERSION).tar.gz fousse@login.medicis.polytechnique.fr:tarball

crq-$(VERSION).tar.gz: $(DIST)
	rm -rf crq-$(VERSION)
	mkdir crq-$(VERSION)
	cp $(DIST) crq-$(VERSION)
	tar czf $@ crq-$(VERSION)
	rm -rf crq-$(VERSION)
