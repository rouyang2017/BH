PROGRAM BH
! Optimization of structure in the coordinate space with BH powered by VASP
! Created by Runhai Ouyang, March, 2023

implicit none
integer nstep,i,j,k,l,natom_tot,ntype,mstep
real*8 prob,rn,dE,temperature,beta,markovE,latt(3,3),energy,mind,maxd,maxdisp(3),tmp
character line*100,command*100,header(9)*100
logical converged,x_dir,y_dir,z_dir
real*8,allocatable:: coord(:,:),oldcoord(:,:)
character,allocatable:: dynamics(:,:)*4,move(:,:)*1

!------------------------------------------------------------
!  Setting your parameters here    
!------------------------------------------------------------
command='srun vasp_std'
natom_tot=85  ! number of atoms in the system
ntype=2       ! number of atom types in POSCAR
nstep=10000  ! total number of MC steps
temperature=3000  ! for Metropolis MC
beta=11604.8/temperature    ! 1/kT
maxdisp=(/0.5,0.5,0.2/)     ! max random atomic displacement (angstrom)
!------------------------------------------------------------

allocate(coord(natom_tot,3))
allocate(oldcoord(natom_tot,3))
allocate(dynamics(natom_tot,3))
allocate(move(natom_tot,3))

mind=0.7      ! interatomic distance has to be > mind*distance
maxd=1.3      ! interatomic distance has to be < maxd*distance
markovE=0.d0   

call random_seed()
open(1,file='BH.out',status='replace')
write(1,'(a)') 'Optimization of structure in the coordinate space with BH powered by VASP'
write(1,'(/a)') 'Input parameters: '
write(1,'(2a)') 'command = ',trim(command)
write(1,'(a,i5)') 'natom_tot = ',natom_tot
write(1,'(a,i5)') 'ntype = ',ntype
write(1,'(a,i10)') 'nstep = ',nstep
write(1,'(a,f10.3)') 'temperature = ',temperature
write(1,'(a,3f10.3)') 'maxdisp = ',maxdisp

! Read the initial POSCAR
open(2,file='POSCAR',status='old')
do i=1,9  
  read(2,'(a)') header(i)
end do
header(9)='Cartesian' 
do i=1,natom_tot
   read(2,'(a)') line
   read(line,*) tmp,tmp,tmp,dynamics(i,:3),move(i,:3)
end do
close(2)


!------------------------------------------
do i=1,nstep

   
   !***************************************
   ! energy calculation, by DFT or by MLFF
   !***************************************
    if(i==1)  call system('mkdir Calcs; cp POSCAR POSCAR_next; cp ML_FF INCAR KPOINTS POTCAR Calcs')
    call system('cd Calcs; mv ../POSCAR_next POSCAR; '//trim(adjustl(command))//'; cd ..')
    call sleep(5) 

   !----------------------------------------------
   ! collect the converged coordinates and energy
   !----------------------------------------------
   converged=.false.
   coord=0.d0
   energy=0.d0
   latt=0.d0

   if(i==1) write(1,'(/a)') ' INCAR: '

   open(2,file='Calcs/OUTCAR',status='old')
   mstep=0
   k=0
   do while(.not. eof(2))
     read(2,'(a)') line
     if(index(line,'direct lattice vectors')/=0) mstep=mstep+1
     if(index(line,'reached required accuracy - stopping structural energy minimisation')/=0) converged=.true.
     if(i==1) then
       if(index(line,'POTCAR:')/=0 .and. k<ntype) then
           write(1,'(a)') trim(line)
           k=k+1
       end if
       if(mstep==2 .and. index(line,'direct lattice vectors')/=0) then
          write(1,'(a)') trim(line)
          do j=1,3
             read(2,'(a)') line
             write(1,'(a)') trim(line)
          end do
       end if
       if(mstep==2 .and. index(line,'position of ions in fractional coordinates')/=0) then
          write(1,'(a)') '  position of ions in fractional coordinates'
          do j=1,natom_tot
             read(2,'(a)') line
             write(1,'(a)') trim(line)
          end do
       end if
       if(index(line,' ions per type')/=0) then
           write(1,'(/a/)') trim(line)
       end if
     end if
   end do

   write(1,'(/a,i6.6/)'), 'Step ',i

   rewind(2)
   k=0
   do while(.not. eof(2))
      read(2,'(a)') line
      if (index(line,'direct lattice vectors')/=0) then
          k=k+1
          if (k==mstep) then
              do l=1,3
                 read(2,*) latt(l,:)
              end do
          end if
      end if
      if(k==mstep .and. index(line,'POSITION')/=0) then
         read(2,*)
         do l=1,natom_tot
            read(2,*) coord(l,:)
         end do
      end if
      if(k==mstep .and. index(line,'energy(sigma->0) =')/=0 ) then
        l=index(line,'energy(sigma->0) =')
        read(line(l+18:),*) energy
      end if
   end do
   close(2)
   
   !*********************************
   ! Metropolis MC selection
   !*********************************
   call random_number(rn) 
   dE=energy-markovE   ! this_energy - previous_energy
   prob=exp(-dE*beta)
   if(prob > rn .and. energy<0.0 ) then
         write(1,'(a,i6.6,a,f20.10,a)') 'Step: ',i,'  Total energy: ',energy,'   Accepted!'
         write (1,'(a,f8.4,a,f8.4)'),'Probability ',prob, ' > random r=',rn
        !-----------------------------------------------------------------
        ! output lattice vectors, coordinates, energy, and convergence
        !-----------------------------------------------------------------
        write(1,'(a)') '  direct lattice vectors'
        do j=1,3
           write(1,'(3f16.9)') latt(j,:)
        end do
        write(1,'(/a)') '  POSITION'
        write(1,'(a)') ' -----------------------------------------------------------------------------------'
        do j=1,natom_tot
           write(1,'(3f13.5)') coord(j,:)
        end do
        write(1,'(a)') ' -----------------------------------------------------------------------------------'
        write(1,'(a,l)') 'Structure converged? ',converged

         markovE=energy
         oldcoord=coord
   else
         write (1,'(a,f8.4,a,f8.4)'),'Probability ',prob, ' < random r=',rn
         write(1,'(a,i6.6,a,f20.10,a)') 'Step: ',i,'  Total energy: ',energy,'   Rejected!'
         coord=oldcoord
   endif

   if(i/=nstep) then
      call  POSCAR_update
   end if

end do

close(1)

deallocate(coord)
deallocate(oldcoord)
deallocate(dynamics)
deallocate(move)

contains

subroutine POSCAR_update
integer i,j,k
real*8 newcoord(natom_tot,3),rand(3)

newcoord=coord
do i=1,natom_tot
   k=0
   do while(k<=10000000)
      k=k+1
      if(k==10000000) print *,'Not sccessful after 10000000 trails!'
      call random_number(rand)
      rand=(rand-0.5)*2
      do j=1,3
          if(index(dynamics(i,j),'T')/=0 .and. move(i,j)=='Y') newcoord(i,j)=coord(i,j)+rand(j)*maxdisp(j)
      end do
      if(isreasonable(newcoord,i)) exit
   end do
end do

! create the new POSCAR with the newcoord
open(10,file='POSCAR_next',status='replace')
do i=1,9
   write(10,'(a)') trim(header(i))
end do

do i=1,natom_tot
   write(10,'(3f20.10,4a)') newcoord(i,:),'    ',dynamics(i,:)
end do

write(1,'(/a)') 'New configuration generated !'
do i=1,natom_tot
   write(1,'(3f20.10)') newcoord(i,:)
end do
write(1,'(a,f10.5)') "Max displacement: ",maxval(abs(coord-newcoord))

close(10)

end subroutine POSCAR_update


function isreasonable(newcoord,n)
logical isreasonable
integer i,j,k,n,natom_typ(ntype)
real*8 newcoord(:,:),dist,minlength,coord_period(9,3)
character element(ntype)*2,atom(natom_tot)*2

isreasonable=.false.
read(header(6),*) element(:ntype)
read(header(7),*) natom_typ(:ntype)

! symbol for each atom
j=1
do i=1,natom_tot
   if(i<=sum(natom_typ(:j))) then
      atom(i)=trim(element(j))
      if(i==sum(natom_typ(:j))) j=j+1
   end if
end do


!Considering periodicity
coord_period(1,:3)=0.d0
coord_period(2,:3)=latt(1,:3)
coord_period(3,:3)=-latt(1,:3)
coord_period(4,:3)=latt(2,:3)
coord_period(5,:3)=-latt(2,:3)
coord_period(6,:3)=latt(1,:3)+latt(2,:3)
coord_period(7,:3)=latt(1,:3)-latt(2,:3)
coord_period(8,:3)=-latt(1,:3)+latt(2,:3)
coord_period(9,:3)=-latt(1,:3)-latt(2,:3)

minlength=10.d0
do k=1,9
   do i=1,natom_tot
      if(i/=n) then
          if(abs(newcoord(i,1)+coord_period(k,1)-newcoord(n,1))>5.0) cycle   !quick screening
          if(abs(newcoord(i,2)+coord_period(k,2)-newcoord(n,2))>5.0) cycle
          if(abs(newcoord(i,3)+coord_period(k,3)-newcoord(n,3))>5.0) cycle

          dist=sqrt(sum((newcoord(i,:)+coord_period(k,:3)-newcoord(n,:))**2))
          dist=dist/radii_sum(atom(i),atom(n))         ! distance normalized by radii_sum
          if(dist<mind) return
          if(dist<minlength) minlength=dist
      end if
   end do
end do

if(minlength>maxd) return
isreasonable=.true.

end function


function radii_sum(atom1,atom2)
! covalent radii(Beatriz Cordero et al., 2018)
integer i,j
character atom1*2,atom2*2,symbol(96)*2
real*8 radii(96),radii_sum

symbol(:96)=(/'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc',&
'Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb',&
'Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd',&
'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au',&
'Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm'/)
radii(:96)=(/0.31,0.28,1.28,0.96,0.84,0.76,0.71,0.66,0.57,0.58,1.66,1.41,1.21,1.11,1.07,1.05,1.02,1.06,2.03,1.76,&
1.70,1.60,1.53,1.39,1.61,1.52,1.50,1.24,1.32,1.22,1.22,1.20,1.19,1.20,1.20,1.16,2.20,1.95,1.90,1.75,1.64,1.54,1.47,&
1.46,1.42,1.39,1.45,1.44,1.42,1.39,1.39,1.38,1.39,1.40,2.44,2.15,2.07,2.04,2.03,2.01,1.99,1.98,1.98,1.96,1.94,1.92,&
1.92,1.89,1.90,1.87,1.87,1.75,1.70,1.62,1.51,1.44,1.41,1.36,1.36,1.32,1.45,1.46,1.48,1.40,1.50,1.50,2.60,2.21,&
2.15,2.06,2.00,1.96,1.90,1.87,1.80,1.69/)

do i=1,96
   do j=1,96
       if(trim(atom1)==trim(symbol(i)) .and. trim(atom2)==trim(symbol(j))) radii_sum=radii(i)+radii(j)
   end do
end do

end function



end
