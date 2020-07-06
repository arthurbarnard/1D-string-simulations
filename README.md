# 1D string simulations
This code simulates timedynamics of generic 1D nanostrings by discretizing them into masses and springs.

##Initial Run of CNT specific code
In order to test this code, first compile `CNT_run.cpp`	to generate an exe file (Dev-C++ is an option). Then run the following command:
```
CNT_run 10 ./initfiles/inputFile ./outfiles/testOut ./initfile/seedTable.txt
```

In this particular format, it expects that there are a series of different inputFiles with an appended number of the form `<input file>_%03d.txt` as well as a seed table with a series of random seeds for the random number generator.

In the above case, we are using `inputFile_010.txt` located in the `./initfiles/` folder and using the $11^{th}$ seed in the seed table, also located in `./initfiles/`
 

