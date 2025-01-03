!------------------------------------------------------------------------------------
! MODULE: Energy
!> @brief Routine for calculating the total energy of the system
!> @author Anders Bergman, Lars Bergqvist, Johan Hellsvik
!> @details The factor fcinv, used in the energy calculation, transforms all the
!> energy contributions to mRy. This is necessary since the parameters needed for the
!> fields are in Teslas.
!> @author
!> Anders Bergman, Johan Hellsvik, Lars Bergqvist, Jonathan Chico
!> @copyright
!! Copyright (C) 2008-2018 UppASD group
!! This file is distributed under the terms of the
!! GNU General Public License.
!! See http://www.gnu.org/copyleft/gpl.txt
!------------------------------------------------------------------------------------
module Energy
   use Parameters
   use Profiling
   use Constants
   use HamiltonianData
   use InputData, only : ham_inp
   use HamiltonianActions
   use LSF, only : totalenergy_LSF
   use DipoleManager, only : dipole_field_calculation,calc_macro_energy

   implicit none

   type ene_t
   ! Separated energy contributions
   real(dblprec), dimension(:), allocatable :: ene_xc
   real(dblprec), dimension(:), allocatable :: ene_dm
   real(dblprec), dimension(:), allocatable :: ene_sa
   real(dblprec), dimension(:), allocatable :: ene_bq
   real(dblprec), dimension(:), allocatable :: ene_ring   
   real(dblprec), dimension(:), allocatable :: ene_pd
   real(dblprec), dimension(:), allocatable :: ene_chir
   real(dblprec), dimension(:), allocatable :: ene_dip
   real(dblprec), dimension(:), allocatable :: ene_ani
   real(dblprec), dimension(:), allocatable :: ene_lsf
   real(dblprec), dimension(:), allocatable :: ene_ext
   real(dblprec), dimension(:), allocatable :: ene_pair
   real(dblprec), dimension(:), allocatable :: ene_bqdm
   !Bq terms
   real(dblprec), dimension(:), allocatable :: ene_bqfull11
   real(dblprec), dimension(:), allocatable :: ene_bqfull21
   real(dblprec), dimension(:), allocatable :: ene_bqfull22
   real(dblprec), dimension(:), allocatable :: ene_bqfull23
   real(dblprec), dimension(:), allocatable :: ene_bqfull31
   real(dblprec), dimension(:), allocatable :: ene_bqfull32
   real(dblprec), dimension(:), allocatable :: ene_bqfull33
   real(dblprec), dimension(:), allocatable :: ene_bqfull34
   real(dblprec), dimension(:), allocatable :: ene_bqfull35
   real(dblprec), dimension(:), allocatable :: ene_bqfull36
   ! Total energy
   real(dblprec), dimension(:), allocatable :: energy
   end type ene_t

   type(ene_t) :: ene

   public

contains

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: calc_energy
   !> @brief Calculates the total energy and the term projected contributions to the energy
   !> @details This routine makes use of the subroutines defined in the effective_field()
   !> this should allow for greater consistency in the calculations.
   !> @author Jonathan Chico
   !---------------------------------------------------------------------------------
   subroutine calc_energy(nHam,mstep,Natom,Nchmax,         &
         conf_num,Mensemble,stop_atom,Num_macro,start_atom, &
      plotenergy,Temp,delta_t,do_lsf,lsf_field,    &
      lsf_interpolate,real_time_measure,simid,cell_index,macro_nlistsize,mmom,emom,  &
      emomM,emomM_macro,external_field,time_external_field,max_no_constellations,    &
      maxNoConstl,unitCellType,constlNCoup,constellations,OPT_flag,                  &
      constellationsNeighType,totene,NA,N1,N2,N3)

      implicit none

      integer, intent(in) :: NA           !< Number of atoms in one cell
      integer, intent(in) :: N1           !< Number of cell repetitions in x direction
      integer, intent(in) :: N2           !< Number of cell repetitions in y direction
      integer, intent(in) :: N3           !< Number of cell repetitions in z direction
      integer, intent(in) :: nHam         !< Number of atoms in Hamiltonian
      integer, intent(in) :: mstep        !< Current simulation step
      integer, intent(in) :: Natom        !< Number of atoms in system
      integer, intent(in) :: Nchmax       !< Number of chemical type
      integer, intent(in) :: conf_num     !< number of configurations for LSF
      integer, intent(in) :: Mensemble    !< Number of ensembles
      integer, intent(in) :: stop_atom    !< Atom to end loop for
      integer, intent(in) :: Num_macro    !< Number of macrocells in the system
      integer, intent(in) :: start_atom   !< Atom to start loop for
      integer, intent(in) :: plotenergy   !< Calculate and plot energy (0/1)
      real(dblprec), intent(in) :: Temp         !< Temperature
      real(dblprec), intent(in) :: delta_t      !< Current time step
      character(len=1), intent(in) :: do_lsf    !< Including LSF energy
      character(len=1), intent(in) :: lsf_field          !< LSF field contribution (Local/Total)
      character(len=1), intent(in) :: lsf_interpolate    !< Interpolate LSF or not
      character(len=1), intent(in) :: real_time_measure  !< Display measurements in real time
      character(len=8), intent(in) :: simid              !< Name of simulation
      integer, dimension(Natom), intent(in) :: cell_index            !< Macrocell index for each atom
      integer, dimension(Num_macro), intent(in) :: macro_nlistsize   !< Number of atoms per macrocell
      real(dblprec), dimension(Natom,Mensemble), intent(in) :: mmom     !< Current magnetic moment
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emom   !< Current unit moment vector
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: emomM  !< Current magnetic moment vector
      real(dblprec), dimension(3,Num_macro,Mensemble), intent(in) :: emomM_macro !< The full vector of the macrocell magnetic moment
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: external_field  !< External magnetic field
      real(dblprec), dimension(3,Natom,Mensemble), intent(in) :: time_external_field !< External time-dependent magnetic field
      ! .. Output Variables
      real(dblprec), intent(out) :: totene !< Total energy
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! +++ Optimization Routines Variables +++ !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer, intent(in) :: max_no_constellations ! The maximum (global) length of the constellation matrix
      integer, dimension(:), intent(in) :: maxNoConstl
      integer, dimension(:,:), intent(in) :: unitCellType ! Array of constellation id and classification (core, boundary, or noise) per atom
      real(dblprec), dimension(:,:,:), intent(in) :: constlNCoup
      real(dblprec), dimension(:,:,:), intent(in) :: constellations
      logical, intent(in) :: OPT_flag
      integer, dimension(:,:,:), intent(in) :: constellationsNeighType
      real(dblprec), dimension(3,max_no_constellations,Mensemble) :: beff1_constellations
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! +++ End Region +++ !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !.. Local scalars
      integer :: ii,kk,inttype
      integer :: i_all, i_stat
      real(dblprec) :: energy_dip
      real(dblprec) :: ene_ext_m, ene_ext_s, fcinv,fc
      real(dblprec) :: exc,edm,ebq,ering,edip,eext,epair,ebqdm,epd,eani,echir, esa
      real(dblprec) :: ebqfull11, ebqfull21
      real(dblprec) :: ebqfull22, ebqfull23
      real(dblprec) :: ebqfull31, ebqfull32, ebqfull33, ebqfull34, ebqfull35, ebqfull36
      real(dblprec) :: energy_m, energy_s, ene_ani_m, ene_ani_s, ene_xc_m, ene_xc_s,ene_lsf_m,ene_lsf_s
      real(dblprec) :: ene_dm_m, ene_dm_s, ene_pd_m, ene_pd_s, ene_bqdm_m, ene_bqdm_s, ene_chir_s, ene_sa_m, ene_sa_s
      real(dblprec) :: ene_bq_m, ene_bq_s, ene_ring_s, ene_ring_m, ene_dip_m, ene_dip_s, ene_pair_m,ene_pair_s, ene_chir_m
      real(dblprec) :: ene_bqfull11_s, ene_bqfull11_m, ene_bqfull21_s, ene_bqfull21_m
      real(dblprec) :: ene_bqfull22_s, ene_bqfull22_m, ene_bqfull23_s, ene_bqfull23_m
      real(dblprec) :: ene_bqfull31_s, ene_bqfull31_m, ene_bqfull32_s, ene_bqfull32_m, ene_bqfull33_s, ene_bqfull33_m, ene_bqfull34_s, ene_bqfull34_m
      real(dblprec) :: ene_bqfull35_s, ene_bqfull35_m, ene_bqfull36_s, ene_bqfull36_m
   
      !.. Local arrays
      character(len=30) :: filn
      real(dblprec), dimension(3) :: beff_xc,beff_dm,beff_pair,beff_pd,beff_bqdm,beff_mdip, beff_chir, beff_sa
      real(dblprec), dimension(3) :: beff_bq,beff_ring,beff_dip,beff_tani,beff_ext,beff_ani,beff_cani
      real(dblprec), dimension(3) :: beff_bqfull11, beff_bqfull21
      real(dblprec), dimension(3) :: beff_bqfull22, beff_bqfull23
      real(dblprec), dimension(3) :: beff_bqfull31, beff_bqfull32, beff_bqfull33, beff_bqfull34, beff_bqfull35, beff_bqfull36
      !      
      real(dblprec), dimension(:,:,:), allocatable :: bfield_dip
      real(dblprec), dimension(:,:,:), allocatable :: site_energy

      !.. Executable statements
      if(OPT_flag) call pre_optimize(Natom,Mensemble,max_no_constellations,maxNoConstl,&
         constellations,constlNCoup,beff1_constellations,constellationsNeighType)

      ! Set the values of all the arrays equal to zero for initialization
      call reset_arrays()

      if(plotenergy==2) then
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! For site-dependent energies, one large array is used. The indices are as follows:
         ! 1: Jij, 2: DM, 3: PseudoDip, 4: BiqDM, 5: BiqH, 6: Dipole, 7: Anisotropy, 8: Zeeman. 9: LSF 10: Chiral, 11:Ring
         ! 12:Bqfull11, 13:Bqfull21, 14:Bqfull22, 15:Bqfull23, 16:Bqfull31, 17:Bqfull32, 18:Bqfull33, 19:Bqfull34, 20:Bqfull35, 21:Bqfull36
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !allocate(site_energy(11,Natom,Mensemble),stat=i_stat)
         allocate(site_energy(21,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(site_energy))*kind(site_energy),'site_energy','calc_energy')
         site_energy=0.0_dblprec
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! If one is considering the dipole-dipole interaction one calls the wrapper
      ! for the calculation of the filed
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (ham_inp%do_dip>0) then
         allocate(bfield_dip(3,Natom,Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(bfield_dip))*kind(bfield_dip),'bfield_dip','calc_energy')
         bfield_dip=0.0_dblprec
         energy_dip=0.0_dblprec
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Wrapper for the calculation of the dipole-dipole interaction field
         ! The field is stored in the bfield array which then is passed to the main loop
         ! This is inefficient for the brute-force methods, but it the best way to ensure
         ! that the FFT approaches can be used in an appropriate way
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         call dipole_field_calculation(NA,N1,N2,N3,Natom,ham_inp%do_dip,Num_macro,          &
            Mensemble,stop_atom,start_atom,cell_index,macro_nlistsize,emomM,        &
            emomM_macro,ham%Qdip,ham%Qdip_macro,energy_dip,bfield_dip)
      endif

      fcinv=mub/mry
      fc = mry/mub
      inttype=1
      if (lsf_interpolate=='L') then
         inttype=0
      elseif(lsf_interpolate=='H') then
         inttype=2
      elseif(lsf_interpolate=='I') then
         inttype=3
      endif

      if (do_lsf=='N') then
         do kk=1, Mensemble
            exc   		= 0.0_dblprec
            edm   		= 0.0_dblprec
            esa   		= 0.0_dblprec
            ebq   		= 0.0_dblprec
            ering 		= 0.0_dblprec
            edip  		= 0.0_dblprec
            eext  		= 0.0_dblprec
            epair 		= 0.0_dblprec
            echir 		= 0.0_dblprec
            ebqdm 		= 0.0_dblprec
            epd   		= 0.0_dblprec
            eani  		= 0.0_dblprec
            !Bq 4-2 interactions
            ebqfull11	= 0.0_dblprec
            ebqfull21	= 0.0_dblprec
            ebqfull22	= 0.0_dblprec
            ebqfull23	= 0.0_dblprec
            ebqfull31	= 0.0_dblprec
            ebqfull32	= 0.0_dblprec
            ebqfull33	= 0.0_dblprec
            ebqfull34	= 0.0_dblprec
            ebqfull35	= 0.0_dblprec
            ebqfull36	= 0.0_dblprec

#if ((! defined  __PATHSCALE__) || (! defined __PGIF90__)) && (!_OPENMP < 201307)
            !$omp parallel do default(shared) schedule(static) &
            !$omp& private(ii,beff_xc,beff_dm,beff_sa,beff_pd,beff_bq,beff_ext,beff_dip,beff_ani,beff_cani,beff_tani,beff_pair,beff_bqdm,beff_mdip) &
            !$omp& reduction(+:exc,edm,epair,epd,ebqdm,ebq,edip,eani,eext,esa,beff_chir)
#endif
            do ii=start_atom, stop_atom

               beff_xc     		= 0.0_dblprec
               beff_dm     		= 0.0_dblprec
               beff_sa     		= 0.0_dblprec
               beff_pd     		= 0.0_dblprec
               beff_bq     		= 0.0_dblprec
               beff_ring   		= 0.0_dblprec              
               beff_ext    		= 0.0_dblprec
               beff_dip    		= 0.0_dblprec
               beff_ani    		= 0.0_dblprec
               beff_cani   		= 0.0_dblprec
               beff_tani   		= 0.0_dblprec
               beff_pair   		= 0.0_dblprec
               beff_bqdm   		= 0.0_dblprec
               beff_mdip   		= 0.0_dblprec
               beff_chir   		= 0.0_dblprec
               !Bq 4-2 interactions
               beff_bqfull11	= 0.0_dblprec
               beff_bqfull21	= 0.0_dblprec
               beff_bqfull22	= 0.0_dblprec
               beff_bqfull23	= 0.0_dblprec
               beff_bqfull31	= 0.0_dblprec
               beff_bqfull32	= 0.0_dblprec
               beff_bqfull33	= 0.0_dblprec
               beff_bqfull34	= 0.0_dblprec
               beff_bqfull35	= 0.0_dblprec
               beff_bqfull36	= 0.0_dblprec

               if(ham_inp%do_jtensor/=1) then
                  ! Heisenberg exchange term
                  if(ham_inp%exc_inter=='N') then
                     call heisenberg_field(ii,kk,beff_xc,Natom,Mensemble,OPT_flag,&
                        beff1_constellations,unitCellType,emomM,max_no_constellations)
                     exc=exc+update_ene(emomM(1:3,ii,kk),beff_xc,0.5_dblprec)
                     if(plotenergy==2) site_energy(1,ii,kk)=update_ene(emomM(1:3,ii,kk),beff_xc,0.5_dblprec)
                  else
                     call heisenberg_rescaling_field(ii,kk,beff_xc,Natom,Mensemble,mmom,emomM)
                     exc=exc+update_ene(emomM(1:3,ii,kk),beff_xc,0.5_dblprec)
                     if(plotenergy==2) site_energy(1,ii,kk)=update_ene(emomM(1:3,ii,kk),beff_xc,0.5_dblprec)
                  endif
                  ! Dzyaloshinskii-Moriya term
                  if(ham_inp%do_dm==1) then
                     call dzyaloshinskii_moriya_field(ii, kk, beff_dm,Natom,Mensemble,emomM)
                     edm=edm+update_ene(emomM(1:3,ii,kk),beff_dm,0.5_dblprec)
                     if(plotenergy==2) site_energy(2,ii,kk)=update_ene(emomM(1:3,ii,kk),beff_dm,0.5_dblprec)
                  endif
                  ! Symmetric anisotropic term
                  if(ham_inp%do_sa==1) then
                     call symmetric_anisotropic_field(ii, kk, beff_sa,Natom,Mensemble,emomM)
                     esa=esa+update_ene(emomM(1:3,ii,kk),beff_sa,0.5_dblprec)
                     if(plotenergy==2) site_energy(2,ii,kk)=update_ene(emomM(1:3,ii,kk),beff_sa,0.5_dblprec)
                  endif
                  beff_pair=beff_xc+beff_dm+beff_sa
               else
                  call tensor_field(ii, kk, beff_pair,Natom,Mensemble,emomM)
                  epair=epair+update_ene(emomM(1:3,ii,kk),beff_pair,0.5_dblprec)
                  if(plotenergy==2) site_energy(2,ii,kk)=update_ene(emomM(1:3,ii,kk),beff_pair,0.5_dblprec)
               end if

               ! Pseudo-Dipolar term
               if(ham_inp%do_pd==1) then
                  call pseudo_dipolar_field(ii, kk, beff_pd,Natom,Mensemble,emomM)
                  epd=epd+update_ene(emomM(1:3,ii,kk),beff_pd,0.5_dblprec)
                  if(plotenergy==2) site_energy(3,ii,kk)=update_ene(emomM(1:3,ii,kk),beff_pd,0.5_dblprec)
               endif

               ! BIQDM term
               if(ham_inp%do_biqdm==1) then
                  call dzyaloshinskii_moriya_bq_field(ii, kk, beff_bqdm,Natom,Mensemble,emomM)
                  ebqdm=ebqdm+update_ene(emomM(1:3,ii,kk),beff_bqdm,0.5_dblprec)
                  if(plotenergy==2) site_energy(4,ii,kk)=update_ene(emomM(1:3,ii,kk),beff_bqdm,0.5_dblprec)
               endif

               ! Biquadratic exchange term
               if(ham_inp%do_bq==1) then
                  call biquadratic_field(ii, kk, beff_bq,Natom,Mensemble,emomM)
                  ebq=ebq+update_ene(emomM(1:3,ii,kk),beff_bq,0.25_dblprec)
                  if(plotenergy==2) site_energy(5,ii,kk)=update_ene(emomM(1:3,ii,kk),beff_bq,0.25_dblprec)
               endif

               ! Four-spin ring exchange term
               if(ham_inp%do_ring==1) then
                  call ring_field(ii, kk, beff_ring,Natom,Mensemble,emomM)
                  ering=ering+update_ene(emomM(1:3,ii,kk),beff_ring,0.25_dblprec)
                  if(plotenergy==2) site_energy(11,ii,kk)=update_ene(emomM(1:3,ii,kk),beff_ring,0.25_dblprec)
               endif

               ! Biquadratic exchange term
               if(ham_inp%do_chir==1) then
                  call chirality_field(ii, kk, beff_chir,Natom,Mensemble,emomM)
                  echir=echir+update_ene(emomM(1:3,ii,kk),beff_chir,0.50_dblprec)
                  if(plotenergy==2) site_energy(10,ii,kk)=update_ene(emomM(1:3,ii,kk),beff_chir,0.5_dblprec)
               endif
               
               !!!!!!
               !Biquadratic 4spin-2site H11 exchange
               if(ham_inp%do_bqfull11==1) then
				  call bqfull11_field(ii, kk, beff_bqfull11,Natom,Mensemble,emomM)
				  ebqfull11=ebqfull11+update_ene(emomM(1:3,ii,kk),beff_bqfull11,0.25_dblprec)
				  if(plotenergy==2) site_energy(12,ii,kk)=update_ene(emomM(1:3,ii,kk),beff_bqfull11,0.5_dblprec)
			   endif
			   ! Biquadratic 4spin-2site H21 interaction
			   if(ham_inp%do_bqfull21==1) then
				  call bqfull21_field(ii, kk, beff_bqfull21,Natom,Mensemble,emomM)
				  ebqfull21=ebqfull21+update_ene(emomM(1:3,ii,kk),beff_bqfull21,0.25_dblprec)
				  if(plotenergy==2) site_energy(13,ii,kk)=update_ene(emomM(1:3,ii,kk),beff_bqfull21,0.5_dblprec)
			   endif
			   !Biquadratic 4spin-2site H22 interaction
               if(ham_inp%do_bqfull22==1) then
				  call bqfull22_field(ii, kk, beff_bqfull22,Natom,Mensemble,emomM)
				  ebqfull22=ebqfull22+update_ene(emomM(1:3,ii,kk),beff_bqfull22,0.25_dblprec)
				  if(plotenergy==2) site_energy(14,ii,kk)=update_ene(emomM(1:3,ii,kk),beff_bqfull22,0.5_dblprec)
			   endif
			   !Biquadratic 4spin-2site H23 interaction
               if(ham_inp%do_bqfull23==1) then
				  call bqfull23_field(ii, kk, beff_bqfull23,Natom,Mensemble,emomM)
				  ebqfull23=ebqfull23+update_ene(emomM(1:3,ii,kk),beff_bqfull23,0.25_dblprec)
				  if(plotenergy==2) site_energy(15,ii,kk)=update_ene(emomM(1:3,ii,kk),beff_bqfull23,0.5_dblprec)
			   endif
			   ! Biquadratic 4spin-2site H31 interaction
			   if(ham_inp%do_bqfull31==1) then
				  call bqfull31_field(ii, kk, beff_bqfull31,Natom,Mensemble,emomM)
				  ebqfull31=ebqfull31+update_ene(emomM(1:3,ii,kk),beff_bqfull31,0.25_dblprec)
				  if(plotenergy==2) site_energy(16,ii,kk)=update_ene(emomM(1:3,ii,kk),beff_bqfull31,0.25_dblprec)
			   endif
			   !Biquadratic 4spin-2site H32 interaction
               if(ham_inp%do_bqfull32==1) then
				  call bqfull32_field(ii, kk, beff_bqfull32,Natom,Mensemble,emomM)
				  ebqfull32=ebqfull32+update_ene(emomM(1:3,ii,kk),beff_bqfull32,0.25_dblprec)
				  if(plotenergy==2) site_energy(17,ii,kk)=update_ene(emomM(1:3,ii,kk),beff_bqfull32,0.25_dblprec)
			   endif
			   !Biquadratic 4spin-2site H33 interaction
               if(ham_inp%do_bqfull33==1) then
				  call bqfull33_field(ii, kk, beff_bqfull33,Natom,Mensemble,emomM)
				  ebqfull33=ebqfull33+update_ene(emomM(1:3,ii,kk),beff_bqfull33,0.25_dblprec)
				  if(plotenergy==2) site_energy(18,ii,kk)=update_ene(emomM(1:3,ii,kk),beff_bqfull33,0.25_dblprec)
			   endif
			   !Biquadratic 4spin-2site H34 interaction
               if(ham_inp%do_bqfull34==1) then
				  call bqfull34_field(ii, kk, beff_bqfull34,Natom,Mensemble,emomM)
				  ebqfull34=ebqfull34+update_ene(emomM(1:3,ii,kk),beff_bqfull34,0.25_dblprec)
				  if(plotenergy==2) site_energy(19,ii,kk)=update_ene(emomM(1:3,ii,kk),beff_bqfull34,0.25_dblprec)
			   endif
			   !Biquadratic 4spin-2site H35 interaction
               if(ham_inp%do_bqfull35==1) then
				  call bqfull35_field(ii, kk, beff_bqfull35,Natom,Mensemble,emomM)
				  ebqfull35=ebqfull35+update_ene(emomM(1:3,ii,kk),beff_bqfull35,0.25_dblprec)
				  if(plotenergy==2) site_energy(20,ii,kk)=update_ene(emomM(1:3,ii,kk),beff_bqfull35,0.25_dblprec)
			   endif
			   !Biquadratic 4spin-2site H36 interaction
               if(ham_inp%do_bqfull36==1) then
				  call bqfull36_field(ii, kk, beff_bqfull36,Natom,Mensemble,emomM)
				  ebqfull36=ebqfull36+update_ene(emomM(1:3,ii,kk),beff_bqfull36,0.25_dblprec)
				  if(plotenergy==2) site_energy(21,ii,kk)=update_ene(emomM(1:3,ii,kk),beff_bqfull36,0.25_dblprec)
			   endif
			   !!!!!!

               ! Dipolar energy contribution
               ! Notice that this makes use of the bfield_dip that is previously calculated
               if (ham_inp%do_dip>0) then
                  ! Site-dependent methods
                  if (ham_inp%do_dip.ne.2) then
                     edip=edip+update_ene(emomM(1:3,ii,kk),bfield_dip(1:3,ii,kk),0.5_dblprec)
                     if(plotenergy==2) site_energy(6,ii,kk)=update_ene(emomM(1:3,ii,kk),bfield_dip(1:3,ii,kk),0.5_dblprec)
                  ! Macrocell method
                  else
                     call calc_macro_energy(ii,kk,bfield_dip(1:3,ii,kk),edip,Natom, &
                        Num_macro,Mensemble,cell_index,emomM_macro,macro_nlistsize)
                  endif
               end if

               if (ham_inp%do_anisotropy==1) then
                  ! Anisotropy
                  if (ham%taniso(ii)==1) then
                     ! Uniaxial anisotropy
                     call uniaxial_anisotropy_field(ii, kk, beff_tani,Natom,Mensemble,ham_inp%mult_axis,emomM)
                     eani=eani+update_ene(emomM(1:3,ii,kk),beff_tani,0.5_dblprec)
                     if(plotenergy==2) site_energy(7,ii,kk)=update_ene(emomM(1:3,ii,kk),beff_tani,0.5_dblprec)
                  elseif (ham%taniso(ii)==2) then
                     ! Cubic anisotropy
                     call cubic_anisotropy_field(ii, kk, beff_tani,Natom,Mensemble,ham_inp%mult_axis,emomM)
                     eani=eani+update_ene(emomM(1:3,ii,kk),beff_tani,0.5_dblprec)
                     if(plotenergy==2) site_energy(7,ii,kk)=update_ene(emomM(1:3,ii,kk),beff_tani,0.5_dblprec)
                  elseif (ham%taniso(ii)==7)then
                     ! Uniaxial and cubic anisotropy
                     call uniaxial_anisotropy_field(ii, kk, beff_ani,Natom,Mensemble,ham_inp%mult_axis,emomM)
                     call cubic_anisotropy_field(ii, kk, beff_cani,Natom,Mensemble,ham_inp%mult_axis,emomM)
                     beff_tani=beff_ani+beff_cani*ham%sb(ii)
                     eani=eani+update_ene(emomM(1:3,ii,kk),beff_tani,0.5_dblprec)
                     if(plotenergy==2) site_energy(7,ii,kk)=update_ene(emomM(1:3,ii,kk),beff_tani,0.5_dblprec)
                  endif
               endif
               ! Contribution of the external field to the energy
               beff_ext=time_external_field(1:3,ii,kk)+external_field(1:3,ii,kk)
               eext=eext+update_ene(emomM(1:3,ii,kk),beff_ext,1.0_dblprec)
               if(plotenergy==2) site_energy(8,ii,kk)=update_ene(emomM(1:3,ii,kk),beff_ext,1.0_dblprec)
            end do
#if ((! defined  __PATHSCALE__) || (! defined __PGIF90__)) && (!_OPENMP < 201307)
         !$omp end parallel do
#endif

            ene%ene_xc(kk)=exc
            ene%ene_dm(kk)=edm
            ene%ene_sa(kk)=esa
            ene%ene_pd(kk)=epd
            ene%ene_bq(kk)=ebq
            ene%ene_ring(kk)=ering            
            ene%ene_chir(kk)=echir
            ene%ene_ext(kk)=eext
            ene%ene_ani(kk)=eani
            ene%ene_dip(kk)=edip
            ene%ene_pair(kk)=epair
            ene%ene_bqdm(kk)=ebqdm
            !Bq 4-2 interactions
            ene%ene_bqfull11(kk)=ebqfull11
            ene%ene_bqfull21(kk)=ebqfull21
            ene%ene_bqfull22(kk)=ebqfull22
            ene%ene_bqfull23(kk)=ebqfull23
            ene%ene_bqfull31(kk)=ebqfull31
            ene%ene_bqfull32(kk)=ebqfull32
            ene%ene_bqfull33(kk)=ebqfull33
            ene%ene_bqfull34(kk)=ebqfull34
            ene%ene_bqfull35(kk)=ebqfull35
            ene%ene_bqfull36(kk)=ebqfull36
         end do
      else
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Calculation of the total LSF energy
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         call totalenergy_LSF(Natom,Nchmax,Mensemble,emom,emomM,mmom, &
            plotenergy,external_field,ene%ene_xc,ene%ene_ani,           &
            ene%ene_ext,ene%ene_lsf,ham_inp%exc_inter,inttype,lsf_field,        &
            site_energy)
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Divide to get energy per atom
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ene%ene_xc(:)=ene%ene_xc(:)/(stop_atom-start_atom+1)
      ene%ene_dm(:)=ene%ene_dm(:)/(stop_atom-start_atom+1)
      ene%ene_sa(:)=ene%ene_sa(:)/(stop_atom-start_atom+1)
      ene%ene_pd(:)=ene%ene_pd(:)/(stop_atom-start_atom+1)
      ene%ene_bq(:)=ene%ene_bq(:)/(stop_atom-start_atom+1)
      ene%ene_ring(:)=ene%ene_ring(:)/(stop_atom-start_atom+1)
      ene%ene_ext(:)=ene%ene_ext(:)/(stop_atom-start_atom+1)
      ene%ene_ani(:)=ene%ene_ani(:)/(stop_atom-start_atom+1)
      ene%ene_dip(:)=ene%ene_dip(:)/(stop_atom-start_atom+1)
      ene%ene_pair(:)=ene%ene_pair(:)/(stop_atom-start_atom+1)
      ene%ene_chir(:)=ene%ene_chir(:)/(stop_atom-start_atom+1)
      ene%ene_bqdm(:)=ene%ene_bqdm(:)/(stop_atom-start_atom+1)
      !4spin-2site interactions
      ene%ene_bqfull11(:)=ene%ene_bqfull11(:)/(stop_atom-start_atom+1)
      ene%ene_bqfull21(:)=ene%ene_bqfull21(:)/(stop_atom-start_atom+1)
      ene%ene_bqfull22(:)=ene%ene_bqfull22(:)/(stop_atom-start_atom+1)
      ene%ene_bqfull23(:)=ene%ene_bqfull23(:)/(stop_atom-start_atom+1)
      ene%ene_bqfull31(:)=ene%ene_bqfull31(:)/(stop_atom-start_atom+1)
      ene%ene_bqfull32(:)=ene%ene_bqfull32(:)/(stop_atom-start_atom+1)
      ene%ene_bqfull33(:)=ene%ene_bqfull33(:)/(stop_atom-start_atom+1)
      ene%ene_bqfull34(:)=ene%ene_bqfull34(:)/(stop_atom-start_atom+1)
      ene%ene_bqfull35(:)=ene%ene_bqfull35(:)/(stop_atom-start_atom+1)
      ene%ene_bqfull36(:)=ene%ene_bqfull36(:)/(stop_atom-start_atom+1)
      !
      ene%ene_lsf(:)=0.50_dblprec*ene%ene_lsf(:)/(stop_atom-start_atom+1)
      if(do_lsf=='Y') ene%ene_xc(:)=0.5_dblprec*ene%ene_xc(:)
      ! Divide the total energy per atom
      if (ham_inp%do_jtensor/=1) then
         ene%energy(:)=ene%ene_xc(:)+ene%ene_dm(:)+ene%ene_pd(:)+ene%ene_bq(:)+     &
            ene%ene_ring(:)+ene%ene_ext(:)+ene%ene_ani(:)+ene%ene_dip(:)+           &
            ene%ene_bqdm(:)+ene%ene_lsf(:)+ene%ene_chir(:)+ene%ene_sa(:) + &
            ene%ene_bqfull11(:)+ene%ene_bqfull21(:)+ene%ene_bqfull22(:)+ene%ene_bqfull23(:) + &
            ene%ene_bqfull31(:)+ene%ene_bqfull32(:)+ene%ene_bqfull33(:)+ene%ene_bqfull34(:) + &
            ene%ene_bqfull35(:)+ene%ene_bqfull36(:)
      else
         ene%energy(:)=ene%ene_pair(:)+ene%ene_pd(:)+ene%ene_bq(:)+ene%ene_ring(:)+ & 
         ene%ene_ext(:)+ene%ene_ani(:)+ene%ene_dip(:)+ene%ene_bqdm(:)+              & 
         ene%ene_lsf(:)+ene%ene_chir(:) + &
         ene%ene_bqfull11(:)+ene%ene_bqfull21(:)+ene%ene_bqfull22(:)+ene%ene_bqfull23(:)+ &
         ene%ene_bqfull31(:)+ene%ene_bqfull32(:)+ene%ene_bqfull33(:)+ene%ene_bqfull34(:)+ &
         ene%ene_bqfull35(:)+ene%ene_bqfull36(:)
      endif

      ! Mean and std.dev. of  energies
      if (ham_inp%do_jtensor/=1) then
         call calculate_mean_and_deviation(ene%ene_xc,Mensemble,ene_xc_m,ene_xc_s,fcinv)
         call calculate_mean_and_deviation(ene%ene_dm,Mensemble,ene_dm_m,ene_dm_s,fcinv)
         call calculate_mean_and_deviation(ene%ene_sa,Mensemble,ene_sa_m,ene_sa_s,fcinv)
      else
         call calculate_mean_and_deviation(ene%ene_pair,Mensemble,ene_pair_m,ene_pair_s,fcinv)
      endif
      call calculate_mean_and_deviation(ene%energy,Mensemble,energy_m,energy_s,fcinv)
      call calculate_mean_and_deviation(ene%ene_pd,Mensemble,ene_pd_m,ene_pd_s,fcinv)
      call calculate_mean_and_deviation(ene%ene_bq,Mensemble,ene_bq_m,ene_bq_s,fcinv)
      call calculate_mean_and_deviation(ene%ene_ring,Mensemble,ene_ring_m,ene_ring_s,fcinv)      
      call calculate_mean_and_deviation(ene%ene_ani,Mensemble,ene_ani_m,ene_ani_s,fcinv)
      call calculate_mean_and_deviation(ene%ene_dip,Mensemble,ene_dip_m,ene_dip_s,fcinv)
      call calculate_mean_and_deviation(ene%ene_ext,Mensemble,ene_ext_m,ene_ext_s,fcinv)
      call calculate_mean_and_deviation(ene%ene_lsf,Mensemble,ene_lsf_m,ene_lsf_s,fcinv)
      call calculate_mean_and_deviation(ene%ene_bqdm,Mensemble,ene_bqdm_m,ene_bqdm_s,fcinv)
      !
      call calculate_mean_and_deviation(ene%ene_bqfull11,Mensemble,ene_bqfull11_m,ene_bqfull11_s,fcinv)
      call calculate_mean_and_deviation(ene%ene_bqfull21,Mensemble,ene_bqfull21_m,ene_bqfull21_s,fcinv)
      call calculate_mean_and_deviation(ene%ene_bqfull22,Mensemble,ene_bqfull22_m,ene_bqfull22_s,fcinv)
      call calculate_mean_and_deviation(ene%ene_bqfull23,Mensemble,ene_bqfull23_m,ene_bqfull23_s,fcinv)
      call calculate_mean_and_deviation(ene%ene_bqfull31,Mensemble,ene_bqfull31_m,ene_bqfull31_s,fcinv)
      call calculate_mean_and_deviation(ene%ene_bqfull32,Mensemble,ene_bqfull32_m,ene_bqfull32_s,fcinv)
      call calculate_mean_and_deviation(ene%ene_bqfull33,Mensemble,ene_bqfull33_m,ene_bqfull33_s,fcinv)
      call calculate_mean_and_deviation(ene%ene_bqfull34,Mensemble,ene_bqfull34_m,ene_bqfull34_s,fcinv)
      call calculate_mean_and_deviation(ene%ene_bqfull35,Mensemble,ene_bqfull35_m,ene_bqfull35_s,fcinv)
      call calculate_mean_and_deviation(ene%ene_bqfull36,Mensemble,ene_bqfull36_m,ene_bqfull36_s,fcinv)
      !
      call calculate_mean_and_deviation(ene%ene_chir,Mensemble,ene_chir_m,ene_chir_s,fcinv)

      ! Rescale energies for other use later (Cv)
      ene%energy=ene%energy*fcinv
      ene%ene_xc=ene%ene_xc*fcinv
      ene%ene_lsf=ene%ene_lsf*fcinv
      totene=energy_m

      ! Print to files
      write (filn,'(''totenergy.'',a,''.out'')') trim(simid)
      open(ofileno, file=filn, position="append")

      if (ham_inp%do_jtensor/=1) then
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Print the total energy
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         if (real_time_measure=='Y') then
            if (mstep-1==0) then
               write(ofileno,10010) "#Time","Tot", "Exc","Ani","DM","PD","BiqDM",   &
               "BQ","Dip","Zeeman","LSF","Chir","Ring","SA","H11(4-2)","H21(4-2)","H22(4-2)","H23(4-2)", &
               "H31(4-2)","H32(4-2)","H33(4-2)","H34(4-2)","H35(4-2)","H36(4-2)"
            endif
            write(ofileno,10005) (mstep-1)*delta_t,energy_m,ene_xc_m,ene_ani_m,     &
               ene_dm_m,ene_pd_m,ene_bqdm_m,ene_bq_m,ene_dip_m,ene_ext_m,ene_lsf_m, &
               ene_chir_m,ene_ring_m,ene_sa_m,&
               ene_bqfull11_m,ene_bqfull21_m,ene_bqfull22_m,ene_bqfull23_m, &
               ene_bqfull31_m,ene_bqfull32_m,ene_bqfull33_m,ene_bqfull34_m, &
               ene_bqfull35_m,ene_bqfull36_m
         else
            if (mstep-1==0) then
               write(ofileno,10010) "#Iter","Tot","Exc","Ani","DM","PD","BiqDM",    &
               "BQ","Dip","Zeeman","LSF","Chir","Ring","SA", "H11(4-2)","H21(4-2)","H22(4-2)","H23(4-2)", &
               "H31(4-2)","H32(4-2)","H33(4-2)","H34(4-2)","H35(4-2)","H36(4-2)"
            endif
            write(ofileno,10004) mstep-1,energy_m,ene_xc_m,ene_ani_m,ene_dm_m,      &
               ene_pd_m,ene_bqdm_m,ene_bq_m,ene_dip_m,ene_ext_m,ene_lsf_m,          &
               ene_chir_m,ene_ring_m,ene_sa_m,&
               ene_bqfull11_m,ene_bqfull21_m,ene_bqfull22_m,ene_bqfull23_m, &
               ene_bqfull31_m,ene_bqfull32_m,ene_bqfull33_m,ene_bqfull34_m, &
               ene_bqfull35_m,ene_bqfull36_m
         endif
         close(ofileno)
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Print the standard deviation of the total energy
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         if (Mensemble>1) then
            write (filn,'(''stdenergy.'',a,''.out'')') trim(simid)
            open(ofileno, file=filn, position="append")
            if (real_time_measure=='Y') then
               if (mstep-1==0) then
                  write(ofileno,10010) "#Time","Tot","Exc","Ani", "DM","PD","BiqDM",&
                  "BQ","Dip Ene","Zeeman","LSF","Chir","Ring", "SA", "H11(4-2)","H21(4-2)","H22(4-2)","H23(4-2)", &
                  "H31(4-2)","H32(4-2)","H33(4-2)","H34(4-2)","H35(4-2)","H36(4-2)"
               endif
               write(ofileno,10005) (mstep-1)*delta_t,energy_s,ene_xc_s,ene_ani_s,  &
                  ene_dm_s,ene_pd_s,ene_bqdm_s,ene_bq_s,ene_dip_s,ene_ext_s,        &
                  ene_lsf_s,ene_chir_s,ene_ring_s, ene_sa_s,&
                  ene_bqfull11_s,ene_bqfull21_s,ene_bqfull22_s,ene_bqfull23_s, &
                  ene_bqfull31_s,ene_bqfull32_s,ene_bqfull33_s,ene_bqfull34_s, &
                  ene_bqfull35_s,ene_bqfull36_s
            else
               if (mstep-1==0) then
                  write(ofileno,10010) "#Iter","Tot","Exc","Ani", "DM","PD","BiqDM",&
                  "BQ","Dip","Zeeman","LSF","Chir","Ring", "SA","H11(4-2)","H21(4-2)","H22(4-2)","H23(4-2)", &
                  "H31(4-2)","H32(4-2)","H33(4-2)","H34(4-2)","H35(4-2)","H36(4-2)"
               endif
               write(ofileno,10004) mstep-1,energy_s,ene_xc_s,ene_ani_s,ene_dm_s,   &
                  ene_pd_s,ene_bqdm_s,ene_bq_s,ene_dip_s,ene_ext_s,ene_lsf_s,       &
                  ene_chir_s,ene_ring_s, ene_sa_s,&
                  ene_bqfull11_s,ene_bqfull21_s,ene_bqfull22_s,ene_bqfull23_s, &
                  ene_bqfull31_s,ene_bqfull32_s,ene_bqfull33_s,ene_bqfull34_s, &
                  ene_bqfull35_s,ene_bqfull36_s
            endif
            close(ofileno)
         endif
      else
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Print the total energy
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         if (real_time_measure=='Y') then
            if (mstep-1==0) then
               write(ofileno,10011) "#Time","Tot","Heis-Tens","Ani","PD","BiqDM",   &
               "BQ","Dip","Zeeman","LSF","Chir","Ring","H11(4-2)","H21(4-2)","H22(4-2)","H23(4-2)", &
               "H31(4-2)","H32(4-2)","H33(4-2)","H34(4-2)","H35(4-2)","H36(4-2)"
            endif
            write(ofileno,10007) (mstep-1)*delta_t,energy_m,ene_pair_m,ene_ani_m,   &
               ene_pd_m,ene_bqdm_m,ene_bq_m,ene_dip_m,ene_ext_m,ene_lsf_m,          &
               ene_chir_m,ene_ring_m,&
               ene_bqfull11_m,ene_bqfull21_m,ene_bqfull22_m,ene_bqfull23_m, &
               ene_bqfull31_m,ene_bqfull32_m,ene_bqfull33_m,ene_bqfull34_m, &
               ene_bqfull35_m,ene_bqfull36_m
         else
            if (mstep-1==0) then
               write(ofileno,10011) "#Iter","Tot","Heis-Tens","Ani","PD","BiqDM",   &
               "BQ","Dip","Zeeman","LSF","Chir","Ring","H11(4-2)","H21(4-2)","H22(4-2)","H23(4-2)", &
               "H31(4-2)","H32(4-2)","H33(4-2)","H34(4-2)","H35(4-2)","H36(4-2)"
            endif
            write(ofileno,10006) mstep-1,energy_m,ene_pair_m,ene_ani_m,ene_pd_m,    &
               ene_bqdm_m,ene_bq_m,ene_dip_m,ene_ext_m,ene_lsf_m,ene_chir_m,        &
               ene_ring_m,&
               ene_bqfull11_m,ene_bqfull21_m,ene_bqfull22_m,ene_bqfull23_m, &
               ene_bqfull31_m,ene_bqfull32_m,ene_bqfull33_m,ene_bqfull34_m, &
               ene_bqfull35_m,ene_bqfull36_m
         endif
         close(ofileno)
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Print the standard deviation of the total energy
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         if (Mensemble>1) then
            write (filn,'(''stdenergy.'',a,''.out'')') trim(simid)
            open(ofileno, file=filn, position="append")
            if (real_time_measure=='Y') then
               if (mstep-1==0) then
                  write(ofileno,10011) "#Time","Tot","Heis-Tens","Ani","PD","BiqDM",&
                  "BQ","Dip","Zeeman","LSF","Chir","Ring","H11(4-2)","H21(4-2)","H22(4-2)","H23(4-2)", &
                  "H31(4-2)","H32(4-2)","H33(4-2)","H34(4-2)","H35(4-2)","H36(4-2)"
               endif
               write(ofileno,10007) (mstep-1)*delta_t,energy_s,ene_pair_s,ene_ani_s,&
               ene_pd_s, ene_bqdm_s,ene_bq_s,ene_dip_s,ene_ext_s,ene_lsf_s,         &
               ene_chir_s,ene_ring_s,&
               ene_bqfull11_s,ene_bqfull21_s,ene_bqfull22_s,ene_bqfull23_s, &
               ene_bqfull31_s,ene_bqfull32_s,ene_bqfull33_s,ene_bqfull34_s, &
               ene_bqfull35_s,ene_bqfull36_s
            else
               if (mstep-1==0) then
                  write(ofileno,10011)"#Iter","Tot","Heis-Tens","Ani","PD","BiqDM",&
                  "BQ","Dipolar","Zeeman","LSF","Chir","Ring","H11(4-2)","H21(4-2)","H22(4-2)","H23(4-2)", &
                  "H31(4-2)","H32(4-2)","H33(4-2)","H34(4-2)","H35(4-2)","H36(4-2)"
               endif
               write(ofileno,10006) mstep-1,energy_s,ene_pair_s,ene_ani_s,ene_pd_s, &
                  ene_bqdm_s,ene_bq_s,ene_dip_s,ene_ext_s,ene_lsf_s,ene_chir_s,     &
                  ene_ring_s,&
                  ene_bqfull11_s,ene_bqfull21_s,ene_bqfull22_s,ene_bqfull23_s, &
                  ene_bqfull31_s,ene_bqfull32_s,ene_bqfull33_s,ene_bqfull34_s, &
                  ene_bqfull35_s,ene_bqfull36_s
            endif
            close(ofileno)
         endif
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Print the local energy contributions
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(plotenergy==2) then
         site_energy=fcinv*site_energy
         write (filn,'(''localenergy.'',a,''.out'')') trim(simid)
         open(ofileno, file=filn, position="append")

         if (mstep-1==0) then
            ! 0: Total, 1: Jij, 2: DM, 3: PseudoDip, 4: BiqDM, 5: BiqH, 6: Dipole, 7: Anisotropy, 8: Zeeman
            ! 12:Bqfull11, 13:Bqfull21, 14:Bqfull22
            write(ofileno,10015) '#Iter','Site','Ens','Tot','Exc','Ani','DM','PD',  &
            'BiqDM','BQ','Dip','Zeeman','LSF','Chir','Ring','H11(4-2)','H21(4-2)','H22(4-2)','H23(4-2)', &
            'H31(4-2)','H32(4-2)','H33(4-2)','H34(4-2)','H35(4-2)','H36(4-2)'
         endif

         do kk=1, Mensemble
            do ii=start_atom, stop_atom
               write(ofileno,10008) mstep-1, ii, kk, sum(site_energy(:,ii,kk)),     &
               site_energy(1,ii,kk),site_energy(7,ii,kk),site_energy(2,ii,kk),      &
               site_energy(3,ii,kk),site_energy(4,ii,kk),site_energy(5,ii,kk),      &
               site_energy(6,ii,kk),site_energy(8,ii,kk),site_energy(9,ii,kk),      &
               site_energy(10,ii,kk),site_energy(11,ii,kk),							&
               site_energy(12,ii,kk),site_energy(13,ii,kk),site_energy(14,ii,kk),site_energy(15,ii,kk), &
               site_energy(16,ii,kk),site_energy(17,ii,kk),site_energy(18,ii,kk),site_energy(19,ii,kk), &
               site_energy(20,ii,kk),site_energy(21,ii,kk)
            end do
         end do
         close(ofileno)

         i_all=-product(shape(site_energy))*kind(site_energy)
         deallocate(site_energy,stat=i_stat)
         call memocc(i_stat,i_all,'site_energy','calc_energy')
      end if

      ! If one considers the dipole-dipole interaction deallocate the respective array
      if (ham_inp%do_dip>0) then
         i_all=-product(shape(bfield_dip))*kind(bfield_dip)
         deallocate(bfield_dip,stat=i_stat)
         call memocc(i_stat,i_all,'bfield_dip','calc_energy')
      endif

      10004 format (i8,24es16.8)
      10005 format (es12.4,24es16.8)
      10006 format (i8,21es16.8)
      10007 format (es12.4,21es16.8)
      10008 format (2i8,i6,22es16.8)
      10010 format (a8,24a20)
      10011 format (a8,21a16)
      10015 format (2a8,a6,22a16)

   end subroutine calc_energy

   !---------------------------------------------------------------------------------
   ! FUNCTION: update_ene
   !> @brief Calculation of the current energy contribution of a given field
   !> @author Jonathan Chico
   !---------------------------------------------------------------------------------
   real(dblprec) function update_ene(curr_emomM,curr_field,factor)

      use Constants

      implicit none

      real(dblprec), intent(in) :: factor
      real(dblprec), dimension(3), intent(in) :: curr_emomM
      real(dblprec), dimension(3), intent(in) :: curr_field

      update_ene=-factor*(sum(curr_emomM(:)*curr_field(:)))

   end function update_ene

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: calculate_mean_and_deviation
   !> @brief Subroutine to calculate the mean and standard deviation of an array
   !---------------------------------------------------------------------------------
   subroutine calculate_mean_and_deviation(vdata,vlen,mean,deviation,factor)
      !
      implicit none
      !
      integer, intent(in) :: vlen
      real(dblprec), dimension(vlen), intent(in) :: vdata
      real(dblprec), intent(out) :: mean
      real(dblprec), intent(out) :: deviation
      real(dblprec), intent(in) :: factor
      !
      if(vlen>1) then
         mean=sum(vdata)
         mean=mean/(vlen*1.0_dblprec)
         deviation=sum((vdata-mean)**2)
         deviation=deviation/((vlen-1.0_dblprec)*1.0_dblprec)
      else
         mean=vdata(1)
         deviation=0.0_dblprec
      end if
      !
      mean=mean*factor
      deviation=deviation*factor
      !
      return
      !
   end subroutine calculate_mean_and_deviation

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: allocate_energies
   !> @brief Subroutine for allocation/deallocation of energy related energy arrays
   !> @author Jonathan Chico
   !---------------------------------------------------------------------------------
   subroutine allocate_energies(flag,Mensemble)

      implicit none

      integer, intent(in) :: flag
      integer, optional, intent(in) :: Mensemble

      integer :: i_stat,i_all

      if (flag>0) then
         allocate(ene%energy(Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(ene%energy))*kind(ene%energy),'ene%energy','allocate_energies')
         ene%energy=0.0_dblprec
         allocate(ene%ene_xc(Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(ene%ene_xc))*kind(ene%ene_xc),'ene%ene_xc','allocate_energies')
         ene%ene_xc=0.0_dblprec
         allocate(ene%ene_dm(Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(ene%ene_dm))*kind(ene%ene_dm),'ene%ene_dm','allocate_energies')
         ene%ene_dm=0.0_dblprec
         allocate(ene%ene_sa(Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(ene%ene_sa))*kind(ene%ene_sa),'ene%ene_sa','allocate_energies')
         ene%ene_sa=0.0_dblprec
         allocate(ene%ene_bq(Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(ene%ene_bq))*kind(ene%ene_bq),'ene%ene_bq','allocate_energies')
         ene%ene_bq=0.0_dblprec
         allocate(ene%ene_ring(Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(ene%ene_ring))*kind(ene%ene_ring),'ene%ene_ring','allocate_energies')
         ene%ene_bq=0.0_dblprec
         allocate(ene%ene_pd(Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(ene%ene_pd))*kind(ene%ene_pd),'ene%ene_pd','allocate_energies')
         ene%ene_pd=0.0_dblprec
         allocate(ene%ene_chir(Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(ene%ene_chir))*kind(ene%ene_chir),'ene%ene_chir','allocate_energies')
         ene%ene_chir=0.0_dblprec
         allocate(ene%ene_ani(Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(ene%ene_ani))*kind(ene%ene_ani),'ene%ene_ani','allocate_energies')
         ene%ene_ani=0.0_dblprec
         allocate(ene%ene_bqdm(Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(ene%ene_bqdm))*kind(ene%ene_bqdm),'ene%ene_bqdm','allocate_energies')
         ene%ene_bqdm=0.0_dblprec
         allocate(ene%ene_ext(Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(ene%ene_ext))*kind(ene%ene_ext),'ene%ene_ext','allocate_energies')
         ene%ene_ext=0.0_dblprec
         allocate(ene%ene_dip(Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(ene%ene_dip))*kind(ene%ene_dip),'ene%ene_dip','allocate_energies')
         ene%ene_dip=0.0_dblprec
         allocate(ene%ene_lsf(Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(ene%ene_lsf))*kind(ene%ene_lsf),'ene%ene_lsf','allocate_energies')
         ene%ene_lsf=0.0_dblprec
         allocate(ene%ene_pair(Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(ene%ene_pair))*kind(ene%ene_pair),'ene%ene_pair','allocate_energies')
         ene%ene_pair=0.0_dblprec
         !!!!!!
         !H11
         allocate(ene%ene_bqfull11(Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(ene%ene_bqfull11))*kind(ene%ene_bqfull11),'ene%ene_bqfull11','allocate_energies')
         ene%ene_bqfull11=0.0_dblprec
         !H21
         allocate(ene%ene_bqfull21(Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(ene%ene_bqfull21))*kind(ene%ene_bqfull21),'ene%ene_bqfull21','allocate_energies')
         ene%ene_bqfull21=0.0_dblprec
         !H22
         allocate(ene%ene_bqfull22(Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(ene%ene_bqfull22))*kind(ene%ene_bqfull22),'ene%ene_bqfull22','allocate_energies')
         ene%ene_bqfull22=0.0_dblprec
         !H23
         allocate(ene%ene_bqfull23(Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(ene%ene_bqfull23))*kind(ene%ene_bqfull23),'ene%ene_bqfull23','allocate_energies')
         ene%ene_bqfull23=0.0_dblprec
         !H31
         allocate(ene%ene_bqfull31(Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(ene%ene_bqfull31))*kind(ene%ene_bqfull31),'ene%ene_bqfull31','allocate_energies')
         ene%ene_bqfull31=0.0_dblprec
         !H32
         allocate(ene%ene_bqfull32(Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(ene%ene_bqfull32))*kind(ene%ene_bqfull32),'ene%ene_bqfull32','allocate_energies')
         ene%ene_bqfull32=0.0_dblprec
         !H33
         allocate(ene%ene_bqfull33(Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(ene%ene_bqfull33))*kind(ene%ene_bqfull33),'ene%ene_bqfull33','allocate_energies')
         ene%ene_bqfull33=0.0_dblprec
         !H34
         allocate(ene%ene_bqfull34(Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(ene%ene_bqfull34))*kind(ene%ene_bqfull34),'ene%ene_bqfull34','allocate_energies')
         ene%ene_bqfull34=0.0_dblprec
         !H35
         allocate(ene%ene_bqfull35(Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(ene%ene_bqfull35))*kind(ene%ene_bqfull35),'ene%ene_bqfull35','allocate_energies')
         ene%ene_bqfull35=0.0_dblprec
         !H36
         allocate(ene%ene_bqfull36(Mensemble),stat=i_stat)
         call memocc(i_stat,product(shape(ene%ene_bqfull36))*kind(ene%ene_bqfull36),'ene%ene_bqfull36','allocate_energies')
         ene%ene_bqfull36=0.0_dblprec
         !!!!!!
      else
         i_all=-product(shape(ene%energy))*kind(ene%energy)
         deallocate(ene%energy,stat=i_stat)
         call memocc(i_stat,i_all,'ene%energy','allocate_energies')
         i_all=-product(shape(ene%ene_xc))*kind(ene%ene_xc)
         deallocate(ene%ene_xc,stat=i_stat)
         call memocc(i_stat,i_all,'ene%ene_xc','allocate_energies')
         i_all=-product(shape(ene%ene_dm))*kind(ene%ene_dm)
         deallocate(ene%ene_dm,stat=i_stat)
         call memocc(i_stat,i_all,'ene%ene_dm','allocate_energies')
         i_all=-product(shape(ene%ene_sa))*kind(ene%ene_sa)
         deallocate(ene%ene_sa,stat=i_stat)
         call memocc(i_stat,i_all,'ene%ene_sa','allocate_energies')
         i_all=-product(shape(ene%ene_bq))*kind(ene%ene_bq)
         deallocate(ene%ene_bq,stat=i_stat)
         call memocc(i_stat,i_all,'ene%ene_bq','allocate_energies')

         i_all=-product(shape(ene%ene_ring))*kind(ene%ene_ring)
         deallocate(ene%ene_ring,stat=i_stat)
         call memocc(i_stat,i_all,'ene%ene_ring','allocate_energies')

         i_all=-product(shape(ene%ene_pd))*kind(ene%ene_pd)
         deallocate(ene%ene_pd,stat=i_stat)
         call memocc(i_stat,i_all,'ene%ene_pd','allocate_energies')
         i_all=-product(shape(ene%ene_chir))*kind(ene%ene_chir)
         deallocate(ene%ene_chir,stat=i_stat)
         call memocc(i_stat,i_all,'ene%ene_chir','allocate_energies')
         i_all=-product(shape(ene%ene_ani))*kind(ene%ene_ani)
         deallocate(ene%ene_ani,stat=i_stat)
         call memocc(i_stat,i_all,'ene%ene_ani','allocate_energies')
         i_all=-product(shape(ene%ene_bqdm))*kind(ene%ene_bqdm)
         deallocate(ene%ene_bqdm,stat=i_stat)
         call memocc(i_stat,i_all,'ene%ene_bqdm','allocate_energies')
         i_all=-product(shape(ene%ene_ext))*kind(ene%ene_ext)
         deallocate(ene%ene_ext,stat=i_stat)
         call memocc(i_stat,i_all,'ene%ene_ext','allocate_energies')
         i_all=-product(shape(ene%ene_dip))*kind(ene%ene_dip)
         deallocate(ene%ene_dip,stat=i_stat)
         call memocc(i_stat,i_all,'ene%ene_dip','allocate_energies')
         i_all=-product(shape(ene%ene_lsf))*kind(ene%ene_lsf)
         deallocate(ene%ene_lsf,stat=i_stat)
         call memocc(i_stat,i_all,'ene%ene_lsf','allocate_energies')
         i_all=-product(shape(ene%ene_pair))*kind(ene%ene_pair)
         deallocate(ene%ene_pair,stat=i_stat)
         call memocc(i_stat,i_all,'ene%ene_pair','allocate_energies')
         !!!!!!
         !H11
         i_all=-product(shape(ene%ene_bqfull11))*kind(ene%ene_bqfull11)
         deallocate(ene%ene_bqfull11,stat=i_stat)
         call memocc(i_stat,i_all,'ene%ene_bqfull11','allocate_energies')
         !H21
         i_all=-product(shape(ene%ene_bqfull21))*kind(ene%ene_bqfull21)
         deallocate(ene%ene_bqfull21,stat=i_stat)
         call memocc(i_stat,i_all,'ene%ene_bqfull21','allocate_energies') 
         !H22
         i_all=-product(shape(ene%ene_bqfull22))*kind(ene%ene_bqfull22)
         deallocate(ene%ene_bqfull22,stat=i_stat)
         call memocc(i_stat,i_all,'ene%ene_bqfull22','allocate_energies') 
         !H23
         i_all=-product(shape(ene%ene_bqfull23))*kind(ene%ene_bqfull23)
         deallocate(ene%ene_bqfull23,stat=i_stat)
         call memocc(i_stat,i_all,'ene%ene_bqfull23','allocate_energies') 
         !H31
         i_all=-product(shape(ene%ene_bqfull31))*kind(ene%ene_bqfull31)
         deallocate(ene%ene_bqfull31,stat=i_stat)
         call memocc(i_stat,i_all,'ene%ene_bqfull31','allocate_energies') 
         !H32
         i_all=-product(shape(ene%ene_bqfull32))*kind(ene%ene_bqfull32)
         deallocate(ene%ene_bqfull32,stat=i_stat)
         call memocc(i_stat,i_all,'ene%ene_bqfull32','allocate_energies') 
         !H33
         i_all=-product(shape(ene%ene_bqfull33))*kind(ene%ene_bqfull33)
         deallocate(ene%ene_bqfull33,stat=i_stat)
         call memocc(i_stat,i_all,'ene%ene_bqfull33','allocate_energies') 
         !H34
         i_all=-product(shape(ene%ene_bqfull34))*kind(ene%ene_bqfull34)
         deallocate(ene%ene_bqfull34,stat=i_stat)
         call memocc(i_stat,i_all,'ene%ene_bqfull34','allocate_energies') 
         !H35
         i_all=-product(shape(ene%ene_bqfull35))*kind(ene%ene_bqfull35)
         deallocate(ene%ene_bqfull35,stat=i_stat)
         call memocc(i_stat,i_all,'ene%ene_bqfull35','allocate_energies')
         !H36
         i_all=-product(shape(ene%ene_bqfull36))*kind(ene%ene_bqfull36)
         deallocate(ene%ene_bqfull36,stat=i_stat)
         call memocc(i_stat,i_all,'ene%ene_bqfull36','allocate_energies') 
         !!!!!!    

      endif

   end subroutine allocate_energies

   !---------------------------------------------------------------------------------
   ! SUBROUTINE: reset_arrays
   !> @brief Subroutine for setting the value of the arrays to be zero for initialization
   !> @author Jonathan Chico
   !---------------------------------------------------------------------------------
   subroutine reset_arrays()

      implicit none

      ene%energy     	= 0.0_dblprec
      ene%ene_xc     	= 0.0_dblprec
      ene%ene_dm     	= 0.0_dblprec
      ene%ene_sa     	= 0.0_dblprec
      ene%ene_pd     	= 0.0_dblprec
      ene%ene_bq     	= 0.0_dblprec
      ene%ene_ring   	= 0.0_dblprec
      ene%ene_chir   	= 0.0_dblprec
      ene%ene_ani    	= 0.0_dblprec
      ene%ene_ext    	= 0.0_dblprec
      ene%ene_dip    	= 0.0_dblprec
      ene%ene_lsf    	= 0.0_dblprec
      ene%ene_pair   	= 0.0_dblprec
      ene%ene_bqdm   	= 0.0_dblprec
      !
      ene%ene_bqfull11	= 0.0_dblprec
      ene%ene_bqfull21	= 0.0_dblprec
      ene%ene_bqfull22	= 0.0_dblprec
      ene%ene_bqfull23	= 0.0_dblprec
      ene%ene_bqfull31	= 0.0_dblprec
      ene%ene_bqfull32	= 0.0_dblprec
      ene%ene_bqfull33	= 0.0_dblprec
      ene%ene_bqfull34	= 0.0_dblprec
      ene%ene_bqfull35	= 0.0_dblprec
      ene%ene_bqfull36	= 0.0_dblprec

   end subroutine reset_arrays

end module Energy
