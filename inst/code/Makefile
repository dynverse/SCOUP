#Makefile
all: scoup scoup_resume cor sp
scoup: scoup.cpp ou.h node.h
	g++ -o scoup scoup.cpp
scoup_resume: scoup_resume.cpp ou.h node.h
	g++ -o scoup_resume scoup_resume.cpp
cor: correlation.cpp ou.h node.h
	g++ -o cor correlation.cpp
sp: sp.cpp sp.h
	g++ -o sp sp.cpp -lblas -llapack
