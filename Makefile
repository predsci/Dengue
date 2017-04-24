FC = gfortran
OBJ = nserovechost.o
FFLAGS = -c -O3 
all: $(OBJ)

#R CMD SHLIB $(OBJ)

nserovechost.o:  nserovechost.f
	$(FC) $(FFLAGS) $(PKG_LIBS) nserovechost.f

#bnldev.o: bnldev.f
#	$(FC) $(FFLAGS) $(PKG_LIBS) bnldev.f

#weibull.o: weibull.f
#	$(FC) $(FFLAGS) $(PKG_LIBS) weibull.f
