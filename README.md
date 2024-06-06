Basin-hopping (BH) program for structure optimization in the coordinate space.   
Ref: https://pubs.acs.org/doi/10.1021/jacs.4c02424

Using the code:
1. Open the BH.f90 and set your parameters at the beginning of the file.
2. Compile the code by: ifort BH.f90 -o BH
3. Put the four VASP files INCAR, POSCAR, POTCAR, KPOINTS under the working directory
4. Put the command ./BH in the submission script and submit the job to run the BH.

Note:
1. Additional information to indicate which atoms to be randomly displaced during BH is needed in the POSCAR (See the provided two examples). The symbol 'Y' denotes Yes (allowing random displacement), and 'N' denotes No.
2. The 'engine' to power the BH can be either DFT or MLFF, depending on the setting in the INCAR.
3. The BH trajectory can be visualized by dragging the output file BH.out into the software Jmol.
