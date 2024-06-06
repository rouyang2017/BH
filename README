Basin-hopping (BH) program for structure optimization in the coordinate space.
Ref: https://pubs.acs.org/doi/10.1021/jacs.4c02424

Using the code:
1. Open the BH.f90 and set your parameters at the beginning of the file.
2. Compile the code: ifort BH.f90 -o BH
3. Put the four VASP files INCAR, POSCAR, POTCAR, KPOINTS under the working directory
4. Put the command ./BH in the submission script and submit the job.

Note:
1. Additional information to indicate which atoms to be randomly displaced is needed in the POSCAR. Two example POSCAR are provided. The symbol 'Y' denotes Yes (allowing random displacement), and 'N' denotes No.
2. The 'engine' to power the BH can be either DFT or MLFF, depending the setting in the INCAR.
3. The output file BH.out can be 
