#MainMakefile.make

RunSimulator:
	make -f CreateSimDirectoriesMakefile.make
	make -f CreateSimDirectoriesMakefile.make run 
	./CreateSimDirectories.sh

clean:
	rm -f CreateSimDirectories.make
	rm -f CreateSimDirectoriesMakefile
	rm -f SimulatePct
