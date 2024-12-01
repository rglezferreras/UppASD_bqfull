!-------------------------------------------------------------------------------
! MODULE: HamiltonianData
!> @brief Data and allocation routines for the Hamiltonian
!> @copyright
!> GNU Public License.
!-------------------------------------------------------------------------------
module HamiltonianData
   use Profiling
   use Parameters
   use HamiltonianDataType
   !
   implicit none
   !
   type(ham_t) :: ham

   public

contains

   !> Allocate arrays for anisotropy
   subroutine allocate_anisotropies(Natom,mult_axis,flag)
      implicit none

      integer, intent(in),optional :: Natom !< Number of atoms in system
      integer, intent(in) :: flag  !< Allocate or deallocate (1/-1)
      character(len=1), intent(in) :: mult_axis !< Flag to treat more than one anisotropy axis at the same time
      integer :: i_all, i_stat

      ! Allocate arrays for anisotropy
      if(flag>0) then
         allocate(ham%taniso(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ham%taniso))*kind(ham%taniso),'taniso','allocate_anisotropies')
         allocate(ham%eaniso(3,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ham%eaniso))*kind(ham%eaniso),'eaniso','allocate_anisotropies')
         allocate(ham%kaniso(2,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ham%kaniso))*kind(ham%kaniso),'kaniso','allocate_anisotropies')
         allocate(ham%sb(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ham%sb))*kind(ham%sb),'sb','allocate_anisotropies')

         if (mult_axis=='Y') then
            allocate(ham%taniso_diff(Natom),stat=i_stat)
            call memocc(i_stat,product(shape(ham%taniso_diff))*kind(ham%taniso_diff),'taniso_diff','allocate_anisotropies')
            allocate(ham%eaniso_diff(3,Natom),stat=i_stat)
            call memocc(i_stat,product(shape(ham%eaniso_diff))*kind(ham%eaniso_diff),'eaniso_diff','allocate_anisotropies')
            allocate(ham%kaniso_diff(2,Natom),stat=i_stat)
            call memocc(i_stat,product(shape(ham%kaniso_diff))*kind(ham%kaniso_diff),'kaniso_diff','allocate_anisotropies')
            allocate(ham%sb_diff(Natom),stat=i_stat)
            call memocc(i_stat,product(shape(ham%sb_diff))*kind(ham%sb_diff),'sb_diff','allocate_anisotropies')
         endif

      else
         i_all=-product(shape(ham%taniso))*kind(ham%taniso)
         deallocate(ham%taniso,stat=i_stat)
         call memocc(i_stat,i_all,'taniso','allocate_anisotropies')
         i_all=-product(shape(ham%eaniso))*kind(ham%eaniso)
         deallocate(ham%eaniso,stat=i_stat)
         call memocc(i_stat,i_all,'eaniso','allocate_anisotropies')
         i_all=-product(shape(ham%kaniso))*kind(ham%kaniso)
         deallocate(ham%kaniso,stat=i_stat)
         call memocc(i_stat,i_all,'kaniso','allocate_anisotropies')
         i_all=-product(shape(ham%sb))*kind(ham%sb)
         deallocate(ham%sb,stat=i_stat)
         call memocc(i_stat,i_all,'sb','allocate_anisotropies')

         if (mult_axis=='Y') then
            i_all=-product(shape(ham%taniso_diff))*kind(ham%taniso_diff)
            deallocate(ham%taniso_diff,stat=i_stat)
            call memocc(i_stat,i_all,'taniso_diff','allocate_anisotropies')
            i_all=-product(shape(ham%eaniso_diff))*kind(ham%eaniso_diff)
            deallocate(ham%eaniso_diff,stat=i_stat)
            call memocc(i_stat,i_all,'eaniso_diff','allocate_anisotropies')
            i_all=-product(shape(ham%kaniso_diff))*kind(ham%kaniso_diff)
            deallocate(ham%kaniso_diff,stat=i_stat)
            call memocc(i_stat,i_all,'kaniso_diff','allocate_anisotropies')
            i_all=-product(shape(ham%sb_diff))*kind(ham%sb_diff)
            deallocate(ham%sb_diff,stat=i_stat)
            call memocc(i_stat,i_all,'sb_diff','allocate_anisotropies')
         endif
      end if

   end subroutine allocate_anisotropies


   !> Allocate arrays for Heisenberg Hamiltonian
   subroutine allocate_hamiltoniandata(Natom,NA,nHam, conf_num,max_no_neigh,do_jtensor,do_lsf,&
         flag,lsf_field,exc_inter)
      implicit none

      integer, intent(in) ::  do_jtensor           !< Use SKKR style exchange tensor (0=off, 1=on)
      integer, intent(in),optional :: NA           !< Number of atoms in unit cell
      integer, intent(in),optional :: nHam         !< Number of atoms in Hamiltonian
      integer, intent(in),optional :: Natom        !< Number of atoms in system
      integer, intent(in),optional :: conf_num     !< Number of configurations for LSF
      integer, intent(in),optional :: max_no_neigh !< Calculated maximum of neighbours for exchange
      character(len=1), intent(in) :: do_lsf       !< Including LSF energy
      character(len=1), intent(in) :: lsf_field    !< LSF field term
      character(len=1), intent(in) :: exc_inter    !< Flag for interpolations of exchange
      integer, intent(in) :: flag  !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat
      ! Exchange
      if(flag>0) then
         allocate(ham%nlistsize(nHam),stat=i_stat)
         call memocc(i_stat,product(shape(ham%nlistsize))*kind(ham%nlistsize),'nlistsize','allocate_hamiltoniandata')
         ham%nlistsize=0
         !ham%NA_red=NA
         !allocate(ham%nlistsize_red(ham%NA_red),stat=i_stat)
         !call memocc(i_stat,product(shape(ham%nlistsize_red))*kind(ham%nlistsize_red),'nlistsize_red','allocate_hamiltoniandata')
         allocate(ham%nlist(max_no_neigh,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ham%nlist))*kind(ham%nlist),'nlist','allocate_hamiltoniandata')
         ham%nlist=0
         allocate(ham%aHam(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ham%aHam))*kind(ham%aHam),'aHam','allocate_hamiltoniandata')
         ham%aHam=0
         if(do_lsf=='Y' .and. lsf_field=='L') then
            allocate(ham%fs_nlistsize(Natom),stat=i_stat)
            call memocc(i_stat,product(shape(ham%fs_nlistsize))*kind(ham%fs_nlistsize),'fs_nlistsize','allocate_hamiltoniandata')
            ham%fs_nlistsize=0
            allocate(ham%fs_nlist(max_no_neigh,Natom),stat=i_stat)
            call memocc(i_stat,product(shape(ham%fs_nlist))*kind(ham%fs_nlist),'fs_nlist','allocate_hamiltoniandata')
            ham%fs_nlist=0
         endif

         if (do_jtensor/=1) then
            allocate(ham%ncoup(max_no_neigh,nHam,conf_num),stat=i_stat)
            call memocc(i_stat,product(shape(ham%ncoup))*kind(ham%ncoup),'ncoup','allocate_hamiltoniandata')
            ham%ncoup=0.0_dblprec
            !allocate(ham%ncoup_red(max_no_neigh,ham%NA_red,conf_num),stat=i_stat)
            !call memocc(i_stat,product(shape(ham%ncoup_red))*kind(ham%ncoup_red),'ncoup_red','allocate_hamiltoniandata')
            if (exc_inter=='Y') then
               allocate(ham%ncoupD(max_no_neigh,nHam,conf_num),stat=i_stat)
               call memocc(i_stat,product(shape(ham%ncoupD))*kind(ham%ncoupD),'ncoupD','allocate_hamiltoniandata')
               ham%ncoupD=0.0_dblprec
               !allocate(ham%ncoupD_red(max_no_neigh,NA,conf_num),stat=i_stat)
               !call memocc(i_stat,product(shape(ham%ncoupD_red))*kind(ham%ncoupD_red),'ncoupD_red','allocate_hamiltoniandata')
            endif
         else
            allocate(ham%j_tens(3,3,max_no_neigh,nHam),stat=i_stat)
            call memocc(i_stat,product(shape(ham%j_tens))*kind(ham%j_tens),'j_tens','allocate_hamiltoniandata')
            ham%j_tens=0.0_dblprec
         end if
      else
         i_all=-product(shape(ham%nlistsize))*kind(ham%nlistsize)
         deallocate(ham%nlistsize,stat=i_stat)
         call memocc(i_stat,i_all,'nlistsize','allocate_hamiltoniandata')
         i_all=-product(shape(ham%nlist))*kind(ham%nlist)
         deallocate(ham%nlist,stat=i_stat)
         call memocc(i_stat,i_all,'nlist','allocate_hamiltoniandata')
         i_all=-product(shape(ham%aHam))*kind(ham%aHam)
         deallocate(ham%aHam,stat=i_stat)
         call memocc(i_stat,i_all,'aHam','allocate_hamiltoniandata')
         if(do_lsf=='Y' .and. lsf_field=='L') then
            i_all=-product(shape(ham%fs_nlistsize))*kind(ham%fs_nlistsize)
            deallocate(ham%fs_nlistsize,stat=i_stat)
            call memocc(i_stat,i_all,'fs_nlistsize','allocate_hamiltoniandata')
            i_all=-product(shape(ham%fs_nlist))*kind(ham%fs_nlist)
            deallocate(ham%fs_nlist,stat=i_stat)
            call memocc(i_stat,i_all,'fs_nlist','allocate_hamiltoniandata')
            i_all=-product(shape(ham%nind))*kind(ham%nind)
            deallocate(ham%nind,stat=i_stat)
            call memocc(i_stat,i_all,'nind','allocate_hamiltoniandata')
         end if
         if (do_jtensor/=1) then
            i_all=-product(shape(ham%ncoup))*kind(ham%ncoup)
            deallocate(ham%ncoup,stat=i_stat)
            call memocc(i_stat,i_all,'ncoup','allocate_hamiltoniandata')
            if (exc_inter=='Y') then
               i_all=-product(shape(ham%ncoupD))*kind(ham%ncoupD)
               deallocate(ham%ncoupD,stat=i_stat)
               call memocc(i_stat,i_all,'ncoupD','allocate_hamiltoniandata')
            endif
         else
            i_all=-product(shape(ham%j_tens))*kind(ham%j_tens)
            deallocate(ham%j_tens,stat=i_stat)
            call memocc(i_stat,i_all,'j_tens','allocate_hamiltoniandata')
         end if
      end if

   end subroutine allocate_hamiltoniandata

   subroutine allocate_hamiltoniandata_ind(flag,Natom,max_no_neigh_ind)

      implicit none

      integer, intent(in) :: flag
      integer, intent(in), optional :: Natom
      integer, intent(in), optional :: max_no_neigh_ind

      integer :: i_stat, i_all

      if (flag>0) then
         allocate(ham%ind_nlistsize(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ham%ind_nlistsize))*kind(ham%ind_nlistsize),'ind_nlistsize','allocate_hamiltoniandata_ind')
         ham%ind_nlistsize=0
         allocate(ham%ind_nlist(max_no_neigh_ind,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ham%ind_nlist))*kind(ham%ind_nlist),'ind_nlist','allocate_hamiltoniandata_ind')
         ham%ind_nlist=0
         allocate(ham%fix_nlistsize(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ham%fix_nlistsize))*kind(ham%fix_nlistsize),'fix_nlistsize','allocate_hamiltoniandata_ind')
         ham%fix_nlistsize=0
         allocate(ham%fix_nlist(max_no_neigh_ind,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ham%fix_nlist))*kind(ham%fix_nlist),'fix_nlist','allocate_hamiltoniandata_ind')
         ham%fix_nlist=0
         allocate(ham%sus_ind(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ham%sus_ind))*kind(ham%sus_ind), 'sus_ind','allocate_hamiltoniandata_ind')
         ham%sus_ind=1.0_dblprec
      else
         i_all=-product(shape(ham%ind_nlistsize))*kind(ham%ind_nlistsize)
         deallocate(ham%ind_nlistsize,stat=i_stat)
         call memocc(i_stat,i_all,'ind_nlistsize','allocate_hamiltoniandata_ind')
         i_all=-product(shape(ham%ind_nlist))*kind(ham%ind_nlist)
         deallocate(ham%ind_nlist,stat=i_stat)
         call memocc(i_stat,i_all,'ind_nlist','allocate_hamiltoniandata_ind')
         i_all=-product(shape(ham%fix_nlistsize))*kind(ham%fix_nlistsize)
         deallocate(ham%fix_nlistsize,stat=i_stat)
         call memocc(i_stat,i_all,'fix_nlistsize','allocate_hamiltoniandata_ind')
         i_all=-product(shape(ham%fix_nlist))*kind(ham%fix_nlist)
         deallocate(ham%fix_nlist,stat=i_stat)
         call memocc(i_stat,i_all,'fix_nlist','allocate_hamiltoniandata_ind')
         i_all=-product(shape(ham%sus_ind))*kind(ham%sus_ind)
         deallocate(ham%sus_ind,stat=i_stat)
         call memocc(i_stat,i_all,'sus_ind','allocate_hamiltoniandata_ind')
         i_all=-product(shape(ham%ind_list_full))*kind(ham%ind_list_full)
         deallocate(ham%ind_list_full,stat=i_stat)
         call memocc(i_stat,i_all,'ind_list_full','allocate_hamiltoniandata_ind')
         if (allocated(ham%fix_list)) then
            i_all=-product(shape(ham%fix_list))*kind(ham%fix_list)
            deallocate(ham%fix_list,stat=i_stat)
            call memocc(i_stat,i_all,'fix_list','allocate_hamiltoniandata_ind')
         endif

      endif

   end subroutine allocate_hamiltoniandata_ind


   !> Allocate arrays for Dzyaloshinskii-Moriya Hamiltonian
   subroutine allocate_dmhamiltoniandata(Natom,nHam, max_no_dmneigh,flag)
      implicit none

      integer, optional, intent(in) :: Natom !< Number of atoms in system
      integer, optional, intent(in) :: nHam !< Number of atoms in Hamiltonian
      integer, optional, intent(in) :: max_no_dmneigh !< Calculated number of neighbours with DM interactions
      integer, intent(in) :: flag !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      if(flag>0) then
         allocate(ham%dmlistsize(nHam),stat=i_stat)
         call memocc(i_stat,product(shape(ham%dmlistsize))*kind(ham%dmlistsize),'dmlistsize','allocate_dmhamiltoniandata')
         allocate(ham%dmlist(max_no_dmneigh,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ham%dmlist))*kind(ham%dmlist),'dmlist','allocate_dmhamiltoniandata')
         allocate(ham%dm_vect(3,max_no_dmneigh,nHam),stat=i_stat)
         call memocc(i_stat,product(shape(ham%dm_vect))*kind(ham%dm_vect),'dm_vect','allocate_dmhamiltoniandata')
      else
         i_all=-product(shape(ham%dmlistsize))*kind(ham%dmlistsize)
         deallocate(ham%dmlistsize,stat=i_stat)
         call memocc(i_stat,i_all,'dmlistsize','allocate_dmhamiltoniandata')
         i_all=-product(shape(ham%dmlist))*kind(ham%dmlist)
         deallocate(ham%dmlist,stat=i_stat)
         call memocc(i_stat,i_all,'dmlist','allocate_dmhamiltoniandata')
         i_all=-product(shape(ham%dm_vect))*kind(ham%dm_vect)
         deallocate(ham%dm_vect,stat=i_stat)
         call memocc(i_stat,i_all,'dm_vect','allocate_dmhamiltoniandata')
      end if

   end subroutine allocate_dmhamiltoniandata

   !> Allocate arrays for Symmetric anisotropic Hamiltonian
   subroutine allocate_sahamiltoniandata(Natom,nHam, max_no_saneigh,flag)
      implicit none

      integer, optional, intent(in) :: Natom !< Number of atoms in system
      integer, optional, intent(in) :: nHam !< Number of atoms in Hamiltonian
      integer, optional, intent(in) :: max_no_saneigh !< Calculated number of neighbours with DM interactions
      integer, intent(in) :: flag !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      if(flag>0) then
         allocate(ham%salistsize(nHam),stat=i_stat)
         call memocc(i_stat,product(shape(ham%salistsize))*kind(ham%salistsize),'salistsize','allocate_sahamiltoniandata')
         allocate(ham%salist(max_no_saneigh,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ham%salist))*kind(ham%salist),'salist','allocate_sahamiltoniandata')
         allocate(ham%sa_vect(3,max_no_saneigh,nHam),stat=i_stat)
         call memocc(i_stat,product(shape(ham%sa_vect))*kind(ham%sa_vect),'sa_vect','allocate_sahamiltoniandata')
      else
         i_all=-product(shape(ham%salistsize))*kind(ham%salistsize)
         deallocate(ham%salistsize,stat=i_stat)
         call memocc(i_stat,i_all,'salistsize','allocate_sahamiltoniandata')
         i_all=-product(shape(ham%salist))*kind(ham%salist)
         deallocate(ham%salist,stat=i_stat)
         call memocc(i_stat,i_all,'salist','allocate_sahamiltoniandata')
         i_all=-product(shape(ham%sa_vect))*kind(ham%sa_vect)
         deallocate(ham%sa_vect,stat=i_stat)
         call memocc(i_stat,i_all,'sa_vect','allocate_sahamiltoniandata')
      end if

   end subroutine allocate_sahamiltoniandata


   !> Allocate arrays for scalar chirality  Hamiltonian
   subroutine allocate_chirhamiltoniandata(Natom,nHam,nn_chir_tot,flag)
      implicit none

      integer, optional, intent(in) :: Natom !< Number of atoms in system
      integer, optional, intent(in) :: nHam  !< Number of atoms in Hamiltonian
      integer, optional, intent(in) :: nn_chir_tot !< Calculated number of neighbours with chir interactions
      integer, intent(in) :: flag !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      if(flag>0) then
         allocate(ham%chirlistsize(nHam),stat=i_stat)
         call memocc(i_stat,product(shape(ham%chirlistsize))*kind(ham%chirlistsize),'chirlistsize','allocate_chirhamiltoniandata')
         allocate(ham%chirlist(2,nn_chir_tot,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ham%chirlist))*kind(ham%chirlist),'chirlist','allocate_chirhamiltoniandata')
         allocate(ham%chir_coup(nn_chir_tot,nHam),stat=i_stat)
         call memocc(i_stat,product(shape(ham%chir_coup))*kind(ham%chir_coup),'chir_coup','allocate_chirhamiltoniandata')
      else
         i_all=-product(shape(ham%chirlistsize))*kind(ham%chirlistsize)
         deallocate(ham%chirlistsize,stat=i_stat)
         call memocc(i_stat,i_all,'chirlistsize','allocate_chirhamiltoniandata')
         i_all=-product(shape(ham%chirlist))*kind(ham%chirlist)
         deallocate(ham%chirlist,stat=i_stat)
         call memocc(i_stat,i_all,'chirlist','allocate_chirhamiltoniandata')
         i_all=-product(shape(ham%chir_coup))*kind(ham%chir_coup)
         deallocate(ham%chir_coup,stat=i_stat)
         call memocc(i_stat,i_all,'chir_coup','allocate_chirhamiltoniandata')
      end if

   end subroutine allocate_chirhamiltoniandata


   !> Allocate arrays for general four-spin Hamiltonian
   subroutine allocate_fourxhamiltoniandata(Natom,nHam,nn_fourx_tot,flag)
      implicit none

      integer, optional, intent(in) :: Natom !< Number of atoms in system
      integer, optional, intent(in) :: nHam  !< Number of atoms in Hamiltonian
      integer, optional, intent(in) :: nn_fourx_tot !< Calculated number of neighbours with fourx interactions
      integer, intent(in) :: flag !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      if(flag>0) then
         allocate(ham%fourxlistsize(nHam),stat=i_stat)
         call memocc(i_stat,product(shape(ham%fourxlistsize))*kind(ham%fourxlistsize),'fourxlistsize','allocate_fourxhamiltoniandata')
         allocate(ham%fourxlist(2,nn_fourx_tot,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ham%fourxlist))*kind(ham%fourxlist),'fourxlist','allocate_fourxhamiltoniandata')
         allocate(ham%fourx_coup(nn_fourx_tot,nHam),stat=i_stat)
         call memocc(i_stat,product(shape(ham%fourx_coup))*kind(ham%fourx_coup),'fourx_coup','allocate_fourxhamiltoniandata')
      else
         i_all=-product(shape(ham%fourxlistsize))*kind(ham%fourxlistsize)
         deallocate(ham%fourxlistsize,stat=i_stat)
         call memocc(i_stat,i_all,'fourxlistsize','allocate_fourxhamiltoniandata')
         i_all=-product(shape(ham%fourxlist))*kind(ham%fourxlist)
         deallocate(ham%fourxlist,stat=i_stat)
         call memocc(i_stat,i_all,'fourxlist','allocate_fourxhamiltoniandata')
         i_all=-product(shape(ham%fourx_coup))*kind(ham%fourx_coup)
         deallocate(ham%fourx_coup,stat=i_stat)
         call memocc(i_stat,i_all,'fourx_coup','allocate_fourxhamiltoniandata')
      end if

   end subroutine allocate_fourxhamiltoniandata


   !> Allocate arrays for Pseudo-Dipolar Hamiltonian
   subroutine allocate_pdhamiltoniandata(Natom,nHam,nn_pd_tot,flag)
      implicit none

      integer, optional, intent(in) :: Natom !< Number of atoms in system
      integer, optional, intent(in) :: nHam  !< Number of atoms in Hamiltonian
      integer, optional, intent(in) :: nn_pd_tot !< Calculated number of neighbours with PD interactions
      integer, intent(in) :: flag !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      if(flag>0) then
         allocate(ham%pdlistsize(nHam),stat=i_stat)
         call memocc(i_stat,product(shape(ham%pdlistsize))*kind(ham%pdlistsize),'pdlistsize','allocate_pdhamiltoniandata')
         allocate(ham%pdlist(nn_pd_tot,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ham%pdlist))*kind(ham%pdlist),'pdlist','allocate_pdhamiltoniandata')
         allocate(ham%pd_vect(9,nn_pd_tot,nHam),stat=i_stat)
         call memocc(i_stat,product(shape(ham%pd_vect))*kind(ham%pd_vect),'pd_vect','allocate_pdhamiltoniandata')
      else
         i_all=-product(shape(ham%pdlistsize))*kind(ham%pdlistsize)
         deallocate(ham%pdlistsize,stat=i_stat)
         call memocc(i_stat,i_all,'pdlistsize','allocate_pdhamiltoniandata')
         i_all=-product(shape(ham%pdlist))*kind(ham%pdlist)
         deallocate(ham%pdlist,stat=i_stat)
         call memocc(i_stat,i_all,'pdlist','allocate_pdhamiltoniandata')
         i_all=-product(shape(ham%pd_vect))*kind(ham%pd_vect)
         deallocate(ham%pd_vect,stat=i_stat)
         call memocc(i_stat,i_all,'pd_vect','allocate_pdhamiltoniandata')
      end if

   end subroutine allocate_pdhamiltoniandata


   !> Allocate arrays for BIQDM Hamiltonian
   subroutine allocate_biqdmhamiltoniandata(Natom,nHam, nn_biqdm_tot,flag)
      implicit none

      integer, optional, intent(in) :: Natom !< Number of atoms in system
      integer, optional, intent(in) :: nHam  !< Number of atoms in Hamiltonian
      integer, optional, intent(in) :: nn_biqdm_tot !< Calculated number of neighbours with BIQDM interactions
      integer, intent(in) :: flag !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      if(flag>0) then
         allocate(ham%biqdmlistsize(nHam),stat=i_stat)
         call memocc(i_stat,product(shape(ham%biqdmlistsize))*kind(ham%biqdmlistsize),'biqdmlistsize','allocate_biqdmhamiltoniandata')
         allocate(ham%biqdmlist(nn_biqdm_tot,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ham%biqdmlist))*kind(ham%biqdmlist),'biqdmlist','allocate_biqdmhamiltoniandata')
         allocate(ham%biqdm_vect(1,nn_biqdm_tot,nHam),stat=i_stat)
         call memocc(i_stat,product(shape(ham%biqdm_vect))*kind(ham%biqdm_vect),'biqdm_vect','allocate_biqdmhamiltoniandata')
      else
         i_all=-product(shape(ham%biqdmlistsize))*kind(ham%biqdmlistsize)
         deallocate(ham%biqdmlistsize,stat=i_stat)
         call memocc(i_stat,i_all,'biqdmlistsize','allocate_biqdmhamiltoniandata')
         i_all=-product(shape(ham%biqdmlist))*kind(ham%biqdmlist)
         deallocate(ham%biqdmlist,stat=i_stat)
         call memocc(i_stat,i_all,'biqdmlist','allocate_biqdmhamiltoniandata')
         i_all=-product(shape(ham%biqdm_vect))*kind(ham%biqdm_vect)
         deallocate(ham%biqdm_vect,stat=i_stat)
         call memocc(i_stat,i_all,'biqdm_vect','allocate_biqdmhamiltoniandata')
      end if

   end subroutine allocate_biqdmhamiltoniandata

   !----------------------------------------------------------------------------
   !> @brief Allocate arrays for biquadratic exchange Hamiltonian
   !----------------------------------------------------------------------------
   subroutine allocate_bqhamiltoniandata(Natom,nHam, nn_bq_tot,flag)
      implicit none

      integer, optional, intent(in) :: Natom !< Number of atoms in system
      integer, optional, intent(in) :: nHam  !< Number of atoms in Hamiltonian
      integer, optional, intent(in) :: nn_bq_tot !< Calculated number of neighbours with BQ interactions
      integer, intent(in) :: flag !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      ! Exchange
      if(flag>0) then
         allocate(ham%bqlistsize(Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ham%bqlistsize))*kind(ham%bqlistsize),'bqlistsize','allocate_bqhamiltoniandata')
         allocate(ham%bqlist(nn_bq_tot,Natom),stat=i_stat)
         call memocc(i_stat,product(shape(ham%bqlist))*kind(ham%bqlist),'bqlist','allocate_bqhamiltoniandata')
         allocate(ham%j_bq(nn_bq_tot,nHam),stat=i_stat)
         call memocc(i_stat,product(shape(ham%j_bq))*kind(ham%j_bq),'j_bq','allocate_bqhamiltoniandata')
      else
         i_all=-product(shape(ham%bqlistsize))*kind(ham%bqlistsize)
         deallocate(ham%bqlistsize,stat=i_stat)
         call memocc(i_stat,i_all,'bqlistsize','allocate_bqhamiltoniandata')
         i_all=-product(shape(ham%bqlist))*kind(ham%bqlist)
         deallocate(ham%bqlist,stat=i_stat)
         call memocc(i_stat,i_all,'bqlist','allocate_bqhamiltoniandata')
         i_all=-product(shape(ham%j_bq))*kind(ham%j_bq)
         deallocate(ham%j_bq,stat=i_stat)
         call memocc(i_stat,i_all,'j_bq','allocate_bqhamiltoniandata')
      end if

   end subroutine allocate_bqhamiltoniandata
   
   !----------------------------------------------------------------------------
   !> @brief Allocate arrays for four-spin ring exchange Hamiltonian
   !----------------------------------------------------------------------------
   
    subroutine allocate_ringhamiltoniandata(Natom,nHam,nn_ring_tot,flag)
   
    implicit none

    integer, optional, intent(in) :: Natom !< Number of atoms in system
    integer, optional, intent(in) :: nHam  !< Number of atoms in Hamiltonian
    integer, optional, intent(in) :: nn_ring_tot !< Calculated number of neighbours with ring interactions
    integer, intent(in) :: flag !< Allocate or deallocate (1/-1)
    integer :: i_all, i_stat

    if(flag>0) then

       allocate(ham%ringlistsize(Natom),stat=i_stat)
       call memocc(i_stat,product(shape(ham%ringlistsize))*kind(ham%ringlistsize),'ringlistsize','allocate_ringhamiltoniandata')
      
       allocate(ham%ringlist(Natom,nn_ring_tot,3),stat=i_stat)
       call memocc(i_stat,product(shape(ham%ringlist))*kind(ham%ringlist),'ringlist','allocate_ringhamiltoniandata')
      
       allocate(ham%j_ring(nHam,nn_ring_tot),stat=i_stat)
       call memocc(i_stat,product(shape(ham%j_ring))*kind(ham%j_ring),'j_ring','allocate_ringhamiltoniandata')

    else

       i_all=-product(shape(ham%ringlistsize))*kind(ham%ringlistsize)
       deallocate(ham%ringlistsize,stat=i_stat)
       call memocc(i_stat,i_all,'ringlistsize','allocate_ringhamiltoniandata')

       i_all=-product(shape(ham%ringlist))*kind(ham%ringlist)
       deallocate(ham%ringlist,stat=i_stat)
       call memocc(i_stat,i_all,'ringlist','allocate_ringhamiltoniandata')

       i_all=-product(shape(ham%j_ring))*kind(ham%j_ring)
       deallocate(ham%j_ring,stat=i_stat)
       call memocc(i_stat,i_all,'j_ring','allocate_ringhamiltoniandata')

    end if

  end subroutine allocate_ringhamiltoniandata

   !----------------------------------------------------------------------------
   !> @brief Allocate arrays for dipole matrix
   !----------------------------------------------------------------------------
   subroutine allocate_dipole(Natom,flag)
      implicit none

      integer, optional, intent(in) :: Natom !< Number of atoms in system
      integer, intent(in) :: flag !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      !  Allocate Q matrix
      if(flag>0) then
         allocate(ham%Qdip(3,3,natom,natom),stat=i_stat)
         call memocc(i_stat,product(shape(ham%Qdip))*kind(ham%Qdip),'Qdip','allocate_dipole')
      else
         i_all=-product(shape(ham%Qdip))*kind(ham%Qdip)
         deallocate(ham%Qdip,stat=i_stat)
         call memocc(i_stat,i_all,'Qdip','allocate_dipole')
      end if

   end subroutine allocate_dipole
   !
   !----------------------------------------------------------------------------
   !> @brief Allocate arrays for macro spin dipole matrix
   !> @author Jonathan Chico
   !----------------------------------------------------------------------------
   subroutine allocate_macro_dipole(Num_macro,flag)

      implicit none

      integer, intent(in) :: Num_macro !< Number of macrocells in the system
      integer, intent(in) :: flag !< Allocate or deallocate (1/-1)

      integer :: i_all, i_stat

      if (flag>0) then
         allocate(ham%Qdip_macro(3,3,Num_macro,Num_macro),stat=i_stat)
         call memocc(i_stat,product(shape(ham%Qdip_macro))*kind(ham%Qdip_macro),'Qdip_macro','allocate_macro_dipole')
         ham%Qdip_macro=0.0_dblprec
      else
         i_all=-product(shape(ham%Qdip_macro))*kind(ham%Qdip_macro)
         deallocate(ham%Qdip_macro,stat=i_stat)
         call memocc(i_stat,i_all,'Qdip_macro','allocate_macro_dipole')
      endif

   end subroutine allocate_macro_dipole


   subroutine scalar_to_tensor(nHam,do_dm,do_sa)
      !
      implicit none
      !
      !
      integer, intent(in) :: nHam  !< Number of atoms in Hamiltonian
      integer, intent(in) :: do_dm   !< Add Dzyaloshinskii-Moriya (DM) term to Hamiltonian (0/1)
      integer, intent(in) :: do_sa   !< Add Symmetric anisotropic (SA) term to Hamiltonian (0/1)
      !
      integer :: i,j, dn, jn, k
      !
      do i=1,nHam
         ham%j_tens(:,:,:,i)=0.0_dblprec
         !Exchange term
         do j=1,ham%nlistsize(i)
            ham%j_tens(1,1,j,i) = ham%ncoup(j,i,1)
            ham%j_tens(2,2,j,i) = ham%ncoup(j,i,1)
            ham%j_tens(3,3,j,i) = ham%ncoup(j,i,1)
         end do

         ! Dzyaloshinskii-Moriya term
         if(do_dm==1) then
            do j=1,ham%dmlistsize(i)
               dn  = ham%dmlist(j,i)
               do k=1,ham%nlistsize(i)
                  jn  = ham%nlist(k,i)
                  if(jn==dn) then
                     ham%j_tens(2,3,j,i) = ham%dm_vect(1,j,i)
                     ham%j_tens(3,2,j,i) = -ham%dm_vect(1,j,i)
                     ham%j_tens(1,3,j,i) = -ham%dm_vect(2,j,i)
                     ham%j_tens(3,1,j,i) = ham%dm_vect(2,j,i)
                     ham%j_tens(1,2,j,i) = ham%dm_vect(3,j,i)
                     ham%j_tens(2,1,j,i) = -ham%dm_vect(3,j,i)
                  end if
               end do
            end do
         end if

         ! Symmetric anisotropic term
         if(do_sa==1) then
            do j=1,ham%salistsize(i)
               dn  = ham%salist(j,i)
               do k=1,ham%nlistsize(i)
                  jn  = ham%nlist(k,i)
                  if(jn==dn) then
                     ham%j_tens(2,3,j,i) = ham%sa_vect(1,j,i)
                     ham%j_tens(3,2,j,i) = ham%sa_vect(1,j,i)
                     ham%j_tens(1,3,j,i) = ham%sa_vect(2,j,i)
                     ham%j_tens(3,1,j,i) = ham%sa_vect(2,j,i)
                     ham%j_tens(1,2,j,i) = ham%sa_vect(3,j,i)
                     ham%j_tens(2,1,j,i) = ham%sa_vect(3,j,i)
                  end if
               end do
            end do
         end if

      end do

   end subroutine scalar_to_tensor
   
   !----------------------------------------------------------------------------
   !> @brief Allocate arrays for biquadratic 4spin-2site interaction H11
   !----------------------------------------------------------------------------
   subroutine allocate_bqfull11data(Natom,nHam,max_no_neigh_bqfull11,flag)
      implicit none

      integer, optional, intent(in) :: Natom 					!< Number of atoms in system
      integer, optional, intent(in) :: nHam 					!< Number of atoms in Hamiltonian
      integer, optional, intent(in) :: max_no_neigh_bqfull11 	! Calculated maximum of neighbours for 4spin-2spin Biquadratic interactions (H11)
      integer, intent(in) 			:: flag						!< Allocate or deallocate (1/-1)
      integer :: i_all, i_stat

      if(flag>0) then
         allocate(ham%bqfull11listsize(nHam),stat=i_stat)
		 call memocc(i_stat,product(shape(ham%bqfull11listsize))*kind(ham%bqfull11listsize),'bqfull11listsize','allocate_bqfull11data')
		 ham%bqfull11listsize=0
		 !
		 allocate(ham%bqfull11list(max_no_neigh_bqfull11,Natom),stat=i_stat)
		 call memocc(i_stat,product(shape(ham%bqfull11list))*kind(ham%bqfull11list),'bqfull11list','allocate_bqfull11data')
		 ham%bqfull11list=0
		 !
		 !allocate(ham%aHam(Natom),stat=i_stat)
		 !call memocc(i_stat,product(shape(ham%aHam))*kind(ham%aHam),'aHam','allocate_hamiltoniandata')
		 !ham%aHam=0
		 !
		 allocate(ham%bqfull11_tens(3,3,max_no_neigh_bqfull11,nHam),stat=i_stat)
		 call memocc(i_stat,product(shape(ham%bqfull11_tens))*kind(ham%bqfull11_tens),'bqfull11_tens','allocate_bqfull11data')
		 ham%bqfull11_tens=0.0_dblprec
      else
		 i_all=-product(shape(ham%bqfull11listsize))*kind(ham%bqfull11listsize)
		 deallocate(ham%bqfull11listsize,stat=i_stat)
		 call memocc(i_stat,i_all,'bqfull11listsize','allocate_bqfull11data')
		 !
		 i_all=-product(shape(ham%bqfull11list))*kind(ham%bqfull11list)
		 deallocate(ham%bqfull11list,stat=i_stat)
		 call memocc(i_stat,i_all,'bqfull11list','allocate_bqfull11data')
		 !
		 !i_all=-product(shape(ham%aHam))*kind(ham%aHam)
		 !deallocate(ham%aHam,stat=i_stat)
		 !call memocc(i_stat,i_all,'aHam','allocate_hamiltoniandata')
		 !
		 i_all=-product(shape(ham%bqfull11_tens))*kind(ham%bqfull11_tens)
		 deallocate(ham%bqfull11_tens,stat=i_stat)
		 call memocc(i_stat,i_all,'bqfull11_tens','allocate_bqfull11data')
      end if
   end subroutine allocate_bqfull11data
   
   
   !----------------------------------------------------------------------------
   !> @brief Allocate arrays for biquadratic 4spin-2site interaction H21
   !----------------------------------------------------------------------------
   subroutine allocate_bqfull21data(Natom,nHam,max_no_neigh_bqfull21,flag)
      implicit none

      integer, optional, intent(in) :: Natom 					!< Number of atoms in system
      integer, optional, intent(in) :: nHam 					!< Number of atoms in Hamiltonian
      integer, optional, intent(in) :: max_no_neigh_bqfull21 	! Calculated maximum of neighbours for 4spin-2spin Biquadratic interactions (H21)
      integer, intent(in) 			:: flag						!< Allocate or deallocate (1/-1)
      integer :: i_all, i_stat

      if(flag>0) then
         allocate(ham%bqfull21listsize(nHam),stat=i_stat)
		 call memocc(i_stat,product(shape(ham%bqfull21listsize))*kind(ham%bqfull21listsize),'bqfull21listsize','allocate_bqfull21data')
		 ham%bqfull21listsize=0
		 !
		 allocate(ham%bqfull21list(max_no_neigh_bqfull21,Natom),stat=i_stat)
		 call memocc(i_stat,product(shape(ham%bqfull21list))*kind(ham%bqfull21list),'bqfull21list','allocate_bqfull21data')
		 ham%bqfull21list=0
		 !
		 !allocate(ham%aHam(Natom),stat=i_stat)
		 !call memocc(i_stat,product(shape(ham%aHam))*kind(ham%aHam),'aHam','allocate_hamiltoniandata')
		 !ham%aHam=0
		 !
		 allocate(ham%bqfull21(max_no_neigh_bqfull21,nHam),stat=i_stat)
		 call memocc(i_stat,product(shape(ham%bqfull21))*kind(ham%bqfull21),'bqfull21','allocate_bqfull21data')
		 ham%bqfull21=0.0_dblprec
      else
		 i_all=-product(shape(ham%bqfull21listsize))*kind(ham%bqfull21listsize)
		 deallocate(ham%bqfull21listsize,stat=i_stat)
		 call memocc(i_stat,i_all,'bqfull21listsize','allocate_bqfull21data')
		 !
		 i_all=-product(shape(ham%bqfull21list))*kind(ham%bqfull21list)
		 deallocate(ham%bqfull21list,stat=i_stat)
		 call memocc(i_stat,i_all,'bqfull21list','allocate_bqfull21data')
		 !
		 !i_all=-product(shape(ham%aHam))*kind(ham%aHam)
		 !deallocate(ham%aHam,stat=i_stat)
		 !call memocc(i_stat,i_all,'aHam','allocate_hamiltoniandata')
		 !
		 i_all=-product(shape(ham%bqfull21))*kind(ham%bqfull21)
		 deallocate(ham%bqfull21,stat=i_stat)
		 call memocc(i_stat,i_all,'bqfull21','allocate_bqfull21data')
      end if
   end subroutine allocate_bqfull21data


   !----------------------------------------------------------------------------
   !> @brief Allocate arrays for biquadratic 4spin-2site interaction H22
   !----------------------------------------------------------------------------
   subroutine allocate_bqfull22data(Natom,nHam,max_no_neigh_bqfull22,flag)
      implicit none

      integer, optional, intent(in) :: Natom 					!< Number of atoms in system
      integer, optional, intent(in) :: nHam 					!< Number of atoms in Hamiltonian
      integer, optional, intent(in) :: max_no_neigh_bqfull22 	! Calculated maximum of neighbours for 4spin-2spin Biquadratic interactions (H22)
      integer, intent(in) 			:: flag 					!< Allocate or deallocate (1/-1)
      integer :: i_all, i_stat
      
      if(flag>0) then
         allocate(ham%bqfull22listsize(nHam),stat=i_stat)
		 call memocc(i_stat,product(shape(ham%bqfull22listsize))*kind(ham%bqfull22listsize),'bqfull22listsize','allocate_bqfull22data')
		 ham%bqfull22listsize=0
		 !
		 allocate(ham%bqfull22list(max_no_neigh_bqfull22,Natom),stat=i_stat)
		 call memocc(i_stat,product(shape(ham%bqfull22list))*kind(ham%bqfull22list),'bqfull22list','allocate_bqfull22data')
		 ham%bqfull22list=0
		 !
		 !allocate(ham%aHam(Natom),stat=i_stat)
		 !call memocc(i_stat,product(shape(ham%aHam))*kind(ham%aHam),'aHam','allocate_hamiltoniandata')
		 !ham%aHam=0
		 !
		 allocate(ham%bqfull22_tens(3,3,max_no_neigh_bqfull22,nHam),stat=i_stat)
		 call memocc(i_stat,product(shape(ham%bqfull22_tens))*kind(ham%bqfull22_tens),'bqfull22_tens','allocate_bqfull22data')
		 ham%bqfull22_tens=0.0_dblprec
      else
		 i_all=-product(shape(ham%bqfull22listsize))*kind(ham%bqfull22listsize)
		 deallocate(ham%bqfull22listsize,stat=i_stat)
		 call memocc(i_stat,i_all,'bqfull22listsize','allocate_bqfull22data')
		 !
		 i_all=-product(shape(ham%bqfull22list))*kind(ham%bqfull22list)
		 deallocate(ham%bqfull22list,stat=i_stat)
		 call memocc(i_stat,i_all,'bqfull22list','allocate_bqfull22data')
		 !
		 !i_all=-product(shape(ham%aHam))*kind(ham%aHam)
		 !deallocate(ham%aHam,stat=i_stat)
		 !call memocc(i_stat,i_all,'aHam','allocate_hamiltoniandata')
		 !
		 i_all=-product(shape(ham%bqfull22_tens))*kind(ham%bqfull22_tens)
		 deallocate(ham%bqfull22_tens,stat=i_stat)
		 call memocc(i_stat,i_all,'bqfull22_tens','allocate_bqfull22data')
      end if
   end subroutine allocate_bqfull22data
   
   !----------------------------------------------------------------------------
   !> @brief Allocate arrays for biquadratic 4spin-2site interaction H23
   !----------------------------------------------------------------------------
   subroutine allocate_bqfull23data(Natom,nHam,max_no_neigh_bqfull23,flag)
      implicit none

      integer, optional, intent(in) :: Natom 					!< Number of atoms in system
      integer, optional, intent(in) :: nHam 					!< Number of atoms in Hamiltonian
      integer, optional, intent(in) :: max_no_neigh_bqfull23 	! Calculated maximum of neighbours for 4spin-2spin Biquadratic interactions (H23)
      integer, intent(in) 			:: flag 					!< Allocate or deallocate (1/-1)
      integer :: i_all, i_stat
    
      if(flag>0) then
         allocate(ham%bqfull23listsize(nHam),stat=i_stat)
		 call memocc(i_stat,product(shape(ham%bqfull23listsize))*kind(ham%bqfull23listsize),'bqfull23listsize','allocate_bqfull23data')
		 ham%bqfull23listsize=0
		 !
		 allocate(ham%bqfull23list(max_no_neigh_bqfull23,Natom),stat=i_stat)
		 call memocc(i_stat,product(shape(ham%bqfull23list))*kind(ham%bqfull23list),'bqfull23list','allocate_bqfull23data')
		 ham%bqfull23list=0
		 !
		 !allocate(ham%aHam(Natom),stat=i_stat)
		 !call memocc(i_stat,product(shape(ham%aHam))*kind(ham%aHam),'aHam','allocate_hamiltoniandata')
		 !ham%aHam=0
		 !
		 allocate(ham%bqfull23_vec(3,max_no_neigh_bqfull23,nHam),stat=i_stat)
		 call memocc(i_stat,product(shape(ham%bqfull23_vec))*kind(ham%bqfull23_vec),'bqfull23_vec','allocate_bqfull23data')
		 ham%bqfull23_vec=0.0_dblprec
      else
		 i_all=-product(shape(ham%bqfull23listsize))*kind(ham%bqfull23listsize)
		 deallocate(ham%bqfull23listsize,stat=i_stat)
		 call memocc(i_stat,i_all,'bqfull23listsize','allocate_bqfull23data')
		 !
		 i_all=-product(shape(ham%bqfull23list))*kind(ham%bqfull23list)
		 deallocate(ham%bqfull23list,stat=i_stat)
		 call memocc(i_stat,i_all,'bqfull23list','allocate_bqfull23data')
		 !
		 !i_all=-product(shape(ham%aHam))*kind(ham%aHam)
		 !deallocate(ham%aHam,stat=i_stat)
		 !call memocc(i_stat,i_all,'aHam','allocate_hamiltoniandata')
		 !
		 i_all=-product(shape(ham%bqfull23_vec))*kind(ham%bqfull23_vec)
		 deallocate(ham%bqfull23_vec,stat=i_stat)
		 call memocc(i_stat,i_all,'bqfull23_vec','allocate_bqfull23data')
      end if
   end subroutine allocate_bqfull23data
   
   
   !----------------------------------------------------------------------------
   !> @brief Allocate arrays for biquadratic 4spin-2site interaction H31
   !----------------------------------------------------------------------------
   subroutine allocate_bqfull31data(Natom,nHam,max_no_neigh_bqfull31,flag)
      implicit none

      integer, optional, intent(in) :: Natom 					!< Number of atoms in system
      integer, optional, intent(in) :: nHam 					!< Number of atoms in Hamiltonian
      integer, optional, intent(in) :: max_no_neigh_bqfull31 	! Calculated maximum of neighbours for 4spin-2spin Biquadratic interactions (H31)
      integer, intent(in) 			:: flag						!< Allocate or deallocate (1/-1)
      integer :: i_all, i_stat

      if(flag>0) then
         allocate(ham%bqfull31listsize(nHam),stat=i_stat)
		 call memocc(i_stat,product(shape(ham%bqfull31listsize))*kind(ham%bqfull31listsize),'bqfull31listsize','allocate_bqfull31data')
		 ham%bqfull31listsize=0
		 !
		 allocate(ham%bqfull31list(max_no_neigh_bqfull31,Natom),stat=i_stat)
		 call memocc(i_stat,product(shape(ham%bqfull31list))*kind(ham%bqfull31list),'bqfull31list','allocate_bqfull31data')
		 ham%bqfull31list=0
		 !
		 !allocate(ham%aHam(Natom),stat=i_stat)
		 !call memocc(i_stat,product(shape(ham%aHam))*kind(ham%aHam),'aHam','allocate_hamiltoniandata')
		 !ham%aHam=0
		 !
		 allocate(ham%bqfull31(max_no_neigh_bqfull31,nHam),stat=i_stat)
		 call memocc(i_stat,product(shape(ham%bqfull31))*kind(ham%bqfull31),'bqfull31','allocate_bqfull31data')
		 ham%bqfull31=0.0_dblprec
      else
		 i_all=-product(shape(ham%bqfull31listsize))*kind(ham%bqfull31listsize)
		 deallocate(ham%bqfull31listsize,stat=i_stat)
		 call memocc(i_stat,i_all,'bqfull31listsize','allocate_bqfull31data')
		 !
		 i_all=-product(shape(ham%bqfull31list))*kind(ham%bqfull31list)
		 deallocate(ham%bqfull31list,stat=i_stat)
		 call memocc(i_stat,i_all,'bqfull31list','allocate_bqfull31data')
		 !
		 !i_all=-product(shape(ham%aHam))*kind(ham%aHam)
		 !deallocate(ham%aHam,stat=i_stat)
		 !call memocc(i_stat,i_all,'aHam','allocate_hamiltoniandata')
		 !
		 i_all=-product(shape(ham%bqfull31))*kind(ham%bqfull31)
		 deallocate(ham%bqfull31,stat=i_stat)
		 call memocc(i_stat,i_all,'bqfull31','allocate_bqfull31data')
      end if
   end subroutine allocate_bqfull31data


   !----------------------------------------------------------------------------
   !> @brief Allocate arrays for biquadratic 4spin-2site interaction H32
   !----------------------------------------------------------------------------
   subroutine allocate_bqfull32data(Natom,nHam,max_no_neigh_bqfull32,flag)
      implicit none

      integer, optional, intent(in) :: Natom 					!< Number of atoms in system
      integer, optional, intent(in) :: nHam 					!< Number of atoms in Hamiltonian
      integer, optional, intent(in) :: max_no_neigh_bqfull32 	! Calculated maximum of neighbours for 4spin-2spin Biquadratic interactions (H32)
      integer, intent(in) 			:: flag 					!< Allocate or deallocate (1/-1)
      integer :: i_all, i_stat
      
      if(flag>0) then
         allocate(ham%bqfull32listsize(nHam),stat=i_stat)
		 call memocc(i_stat,product(shape(ham%bqfull32listsize))*kind(ham%bqfull32listsize),'bqfull32listsize','allocate_bqfull32data')
		 ham%bqfull32listsize=0
		 !
		 allocate(ham%bqfull32list(max_no_neigh_bqfull32,Natom),stat=i_stat)
		 call memocc(i_stat,product(shape(ham%bqfull32list))*kind(ham%bqfull32list),'bqfull32list','allocate_bqfull32data')
		 ham%bqfull32list=0
		 !
		 !allocate(ham%aHam(Natom),stat=i_stat)
		 !call memocc(i_stat,product(shape(ham%aHam))*kind(ham%aHam),'aHam','allocate_hamiltoniandata')
		 !ham%aHam=0
		 !
		 allocate(ham%bqfull32_tens(3,3,max_no_neigh_bqfull32,nHam),stat=i_stat)
		 call memocc(i_stat,product(shape(ham%bqfull32_tens))*kind(ham%bqfull32_tens),'bqfull32_tens','allocate_bqfull32data')
		 ham%bqfull32_tens=0.0_dblprec
      else
		 i_all=-product(shape(ham%bqfull32listsize))*kind(ham%bqfull32listsize)
		 deallocate(ham%bqfull32listsize,stat=i_stat)
		 call memocc(i_stat,i_all,'bqfull32listsize','allocate_bqfull32data')
		 !
		 i_all=-product(shape(ham%bqfull32list))*kind(ham%bqfull32list)
		 deallocate(ham%bqfull32list,stat=i_stat)
		 call memocc(i_stat,i_all,'bqfull32list','allocate_bqfull32data')
		 !
		 !i_all=-product(shape(ham%aHam))*kind(ham%aHam)
		 !deallocate(ham%aHam,stat=i_stat)
		 !call memocc(i_stat,i_all,'aHam','allocate_hamiltoniandata')
		 !
		 i_all=-product(shape(ham%bqfull32_tens))*kind(ham%bqfull32_tens)
		 deallocate(ham%bqfull32_tens,stat=i_stat)
		 call memocc(i_stat,i_all,'bqfull32_tens','allocate_bqfull32data')
      end if
   end subroutine allocate_bqfull32data
   

   !----------------------------------------------------------------------------
   !> @brief Allocate arrays for biquadratic 4spin-2site interaction H33
   !----------------------------------------------------------------------------
   subroutine allocate_bqfull33data(Natom,nHam,max_no_neigh_bqfull33,flag)
      implicit none

      integer, optional, intent(in) :: Natom 					!< Number of atoms in system
      integer, optional, intent(in) :: nHam 					!< Number of atoms in Hamiltonian
      integer, optional, intent(in) :: max_no_neigh_bqfull33 	! Calculated maximum of neighbours for 4spin-2spin Biquadratic interactions (H33)
      integer, intent(in) 			:: flag 					!< Allocate or deallocate (1/-1)
      integer :: i_all, i_stat
      
      if(flag>0) then
         allocate(ham%bqfull33listsize(nHam),stat=i_stat)
		 call memocc(i_stat,product(shape(ham%bqfull33listsize))*kind(ham%bqfull33listsize),'bqfull33listsize','allocate_bqfull33data')
		 ham%bqfull33listsize=0
		 !
		 allocate(ham%bqfull33list(max_no_neigh_bqfull33,Natom),stat=i_stat)
		 call memocc(i_stat,product(shape(ham%bqfull33list))*kind(ham%bqfull33list),'bqfull33list','allocate_bqfull33data')
		 ham%bqfull33list=0
		 !
		 !allocate(ham%aHam(Natom),stat=i_stat)
		 !call memocc(i_stat,product(shape(ham%aHam))*kind(ham%aHam),'aHam','allocate_hamiltoniandata')
		 !ham%aHam=0
		 !
		 allocate(ham%bqfull33_tens(3,3,max_no_neigh_bqfull33,nHam),stat=i_stat)
		 call memocc(i_stat,product(shape(ham%bqfull33_tens))*kind(ham%bqfull33_tens),'bqfull33_tens','allocate_bqfull33data')
		 ham%bqfull33_tens=0.0_dblprec
      else
		 i_all=-product(shape(ham%bqfull33listsize))*kind(ham%bqfull33listsize)
		 deallocate(ham%bqfull33listsize,stat=i_stat)
		 call memocc(i_stat,i_all,'bqfull33listsize','allocate_bqfull33data')
		 !
		 i_all=-product(shape(ham%bqfull33list))*kind(ham%bqfull33list)
		 deallocate(ham%bqfull33list,stat=i_stat)
		 call memocc(i_stat,i_all,'bqfull33list','allocate_bqfull33data')
		 !
		 !i_all=-product(shape(ham%aHam))*kind(ham%aHam)
		 !deallocate(ham%aHam,stat=i_stat)
		 !call memocc(i_stat,i_all,'aHam','allocate_hamiltoniandata')
		 !
		 i_all=-product(shape(ham%bqfull33_tens))*kind(ham%bqfull33_tens)
		 deallocate(ham%bqfull33_tens,stat=i_stat)
		 call memocc(i_stat,i_all,'bqfull33_tens','allocate_bqfull33data')
      end if
   end subroutine allocate_bqfull33data

   !----------------------------------------------------------------------------
   !> @brief Allocate arrays for biquadratic 4spin-2site interaction H34
   !----------------------------------------------------------------------------
   subroutine allocate_bqfull34data(Natom,nHam,max_no_neigh_bqfull34,flag)
      implicit none

      integer, optional, intent(in) :: Natom 					!< Number of atoms in system
      integer, optional, intent(in) :: nHam 					!< Number of atoms in Hamiltonian
      integer, optional, intent(in) :: max_no_neigh_bqfull34 	! Calculated maximum of neighbours for 4spin-2spin Biquadratic interactions (H34)
      integer, intent(in) 			:: flag 					!< Allocate or deallocate (1/-1)
      integer :: i_all, i_stat
    
      if(flag>0) then
         allocate(ham%bqfull34listsize(nHam),stat=i_stat)
		 call memocc(i_stat,product(shape(ham%bqfull34listsize))*kind(ham%bqfull34listsize),'bqfull34listsize','allocate_bqfull34data')
		 ham%bqfull34listsize=0
		 !
		 allocate(ham%bqfull34list(max_no_neigh_bqfull34,Natom),stat=i_stat)
		 call memocc(i_stat,product(shape(ham%bqfull34list))*kind(ham%bqfull34list),'bqfull34list','allocate_bqfull34data')
		 ham%bqfull34list=0
		 !
		 !allocate(ham%aHam(Natom),stat=i_stat)
		 !call memocc(i_stat,product(shape(ham%aHam))*kind(ham%aHam),'aHam','allocate_hamiltoniandata')
		 !ham%aHam=0
		 !
		 allocate(ham%bqfull34_vec(3,max_no_neigh_bqfull34,nHam),stat=i_stat)
		 call memocc(i_stat,product(shape(ham%bqfull34_vec))*kind(ham%bqfull34_vec),'bqfull34_vec','allocate_bqfull34data')
		 ham%bqfull34_vec=0.0_dblprec
      else
		 i_all=-product(shape(ham%bqfull34listsize))*kind(ham%bqfull34listsize)
		 deallocate(ham%bqfull34listsize,stat=i_stat)
		 call memocc(i_stat,i_all,'bqfull34listsize','allocate_bqfull34data')
		 !
		 i_all=-product(shape(ham%bqfull34list))*kind(ham%bqfull34list)
		 deallocate(ham%bqfull34list,stat=i_stat)
		 call memocc(i_stat,i_all,'bqfull34list','allocate_bqfull34data')
		 !
		 !i_all=-product(shape(ham%aHam))*kind(ham%aHam)
		 !deallocate(ham%aHam,stat=i_stat)
		 !call memocc(i_stat,i_all,'aHam','allocate_hamiltoniandata')
		 !
		 i_all=-product(shape(ham%bqfull34_vec))*kind(ham%bqfull34_vec)
		 deallocate(ham%bqfull34_vec,stat=i_stat)
		 call memocc(i_stat,i_all,'bqfull34_vec','allocate_bqfull34data')
      end if
   end subroutine allocate_bqfull34data
 

   !----------------------------------------------------------------------------
   !> @brief Allocate arrays for biquadratic 4spin-2site interaction H35
   !----------------------------------------------------------------------------
   subroutine allocate_bqfull35data(Natom,nHam,max_no_neigh_bqfull35,flag)
      implicit none

      integer, optional, intent(in) :: Natom 					!< Number of atoms in system
      integer, optional, intent(in) :: nHam 					!< Number of atoms in Hamiltonian
      integer, optional, intent(in) :: max_no_neigh_bqfull35 	! Calculated maximum of neighbours for 4spin-2spin Biquadratic interactions (H35)
      integer, intent(in) 			:: flag 					!< Allocate or deallocate (1/-1)
      integer :: i_all, i_stat
    
      if(flag>0) then
         allocate(ham%bqfull35listsize(nHam),stat=i_stat)
		 call memocc(i_stat,product(shape(ham%bqfull35listsize))*kind(ham%bqfull35listsize),'bqfull35listsize','allocate_bqfull35data')
		 ham%bqfull35listsize=0
		 !
		 allocate(ham%bqfull35list(max_no_neigh_bqfull35,Natom),stat=i_stat)
		 call memocc(i_stat,product(shape(ham%bqfull35list))*kind(ham%bqfull35list),'bqfull35list','allocate_bqfull35data')
		 ham%bqfull35list=0
		 !
		 !allocate(ham%aHam(Natom),stat=i_stat)
		 !call memocc(i_stat,product(shape(ham%aHam))*kind(ham%aHam),'aHam','allocate_hamiltoniandata')
		 !ham%aHam=0
		 !
		 allocate(ham%bqfull35_3tens(3,3,3,max_no_neigh_bqfull35,nHam),stat=i_stat)
		 call memocc(i_stat,product(shape(ham%bqfull35_3tens))*kind(ham%bqfull35_3tens),'bqfull35_3tens','allocate_bqfull35data')
		 ham%bqfull35_3tens=0.0_dblprec
      else
		 i_all=-product(shape(ham%bqfull35listsize))*kind(ham%bqfull35listsize)
		 deallocate(ham%bqfull35listsize,stat=i_stat)
		 call memocc(i_stat,i_all,'bqfull35listsize','allocate_bqfull35data')
		 !
		 i_all=-product(shape(ham%bqfull35list))*kind(ham%bqfull35list)
		 deallocate(ham%bqfull35list,stat=i_stat)
		 call memocc(i_stat,i_all,'bqfull35list','allocate_bqfull35data')
		 !
		 !i_all=-product(shape(ham%aHam))*kind(ham%aHam)
		 !deallocate(ham%aHam,stat=i_stat)
		 !call memocc(i_stat,i_all,'aHam','allocate_hamiltoniandata')
		 !
		 i_all=-product(shape(ham%bqfull35_3tens))*kind(ham%bqfull35_3tens)
		 deallocate(ham%bqfull35_3tens,stat=i_stat)
		 call memocc(i_stat,i_all,'bqfull35_3tens','allocate_bqfull35data')
      end if
   end subroutine allocate_bqfull35data

   !----------------------------------------------------------------------------
   !> @brief Allocate arrays for biquadratic 4spin-2site interaction H36
   !----------------------------------------------------------------------------
   subroutine allocate_bqfull36data(Natom,nHam,max_no_neigh_bqfull36,flag)
      implicit none

      integer, optional, intent(in) :: Natom 					!< Number of atoms in system
      integer, optional, intent(in) :: nHam 					!< Number of atoms in Hamiltonian
      integer, optional, intent(in) :: max_no_neigh_bqfull36 	! Calculated maximum of neighbours for 4spin-2spin Biquadratic interactions (H36)
      integer, intent(in) 			:: flag 					!< Allocate or deallocate (1/-1)
      integer :: i_all, i_stat
    
      if(flag>0) then
         allocate(ham%bqfull36listsize(nHam),stat=i_stat)
		 call memocc(i_stat,product(shape(ham%bqfull36listsize))*kind(ham%bqfull36listsize),'bqfull36listsize','allocate_bqfull36data')
		 ham%bqfull36listsize=0
		 !
		 allocate(ham%bqfull36list(max_no_neigh_bqfull36,Natom),stat=i_stat)
		 call memocc(i_stat,product(shape(ham%bqfull36list))*kind(ham%bqfull36list),'bqfull36list','allocate_bqfull36data')
		 ham%bqfull36list=0
		 !
		 !allocate(ham%aHam(Natom),stat=i_stat)
		 !call memocc(i_stat,product(shape(ham%aHam))*kind(ham%aHam),'aHam','allocate_hamiltoniandata')
		 !ham%aHam=0
		 !
		 allocate(ham%bqfull36_3tens(3,3,3,max_no_neigh_bqfull36,nHam),stat=i_stat)
		 call memocc(i_stat,product(shape(ham%bqfull36_3tens))*kind(ham%bqfull36_3tens),'bqfull36_3tens','allocate_bqfull36data')
		 ham%bqfull36_3tens=0.0_dblprec
      else
		 i_all=-product(shape(ham%bqfull36listsize))*kind(ham%bqfull36listsize)
		 deallocate(ham%bqfull36listsize,stat=i_stat)
		 call memocc(i_stat,i_all,'bqfull36listsize','allocate_bqfull36data')
		 !
		 i_all=-product(shape(ham%bqfull36list))*kind(ham%bqfull36list)
		 deallocate(ham%bqfull36list,stat=i_stat)
		 call memocc(i_stat,i_all,'bqfull36list','allocate_bqfull36data')
		 !
		 !i_all=-product(shape(ham%aHam))*kind(ham%aHam)
		 !deallocate(ham%aHam,stat=i_stat)
		 !call memocc(i_stat,i_all,'aHam','allocate_hamiltoniandata')
		 !
		 i_all=-product(shape(ham%bqfull36_3tens))*kind(ham%bqfull36_3tens)
		 deallocate(ham%bqfull36_3tens,stat=i_stat)
		 call memocc(i_stat,i_all,'bqfull36_3tens','allocate_bqfull36data')
      end if
   end subroutine allocate_bqfull36data
end module HamiltonianData
