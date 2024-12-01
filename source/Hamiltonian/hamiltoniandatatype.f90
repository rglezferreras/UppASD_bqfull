!-------------------------------------------------------------------------------
! MODULE: HamiltonianData
!> @brief Data and allocation routines for the Hamiltonian
!> @copyright
!> GNU Public License.
!-------------------------------------------------------------------------------
module HamiltonianDataType

   use Parameters
   !
   implicit none
   !
   ! From setup
   type ham_t
      integer, dimension(:), allocatable :: aHam !< Lookup-table for Hamiltonian
      ! Variables for Heisenberg exchange
      integer ::  max_no_neigh                                 !< Calculated maximum of neighbours for exchange
      integer, dimension(:), allocatable :: nlistsize          !< Size of neighbour list for Heisenberg exchange couplings
      integer, dimension(:,:), allocatable :: nlist            !< Neighbour list for Heisenberg exchange couplings
      real(dblprec), dimension(:,:,:), allocatable :: ncoup    !< Heisenberg exchange couplings
      real(dblprec), dimension(:,:,:), allocatable :: ncoupD   !< Heisenberg exchange couplings (DLM)
      real(dblprec), dimension(:,:,:,:), allocatable :: j_tens !< Exchange tensor (SKKR style)
      !! Reduced variables for Heisenberg exchange
      !integer :: NA_red                                          !< Number of reduced atoms
      !integer, dimension(:), allocatable :: nlistsize_red        !< Size of neighbour list for Heisenberg exchange couplings (reduced)
      !integer, dimension(:,:), allocatable :: nlist_red          !< Neighbour list for Heisenberg exchange couplings (reduced)
      !real(dblprec), dimension(:,:,:), allocatable :: ncoup_red  !< Heisenberg exchange couplings (reduced)
      !real(dblprec), dimension(:,:,:), allocatable :: ncoupD_red !< Heisenberg exchange couplings (DLM) (reduced)
      ! Variables for DMI
      integer ::  max_no_dmneigh                               !< Calculated maximum of neighbours for DM exchange
      integer, dimension(:), allocatable :: dmlistsize         !< Size of neighbour list for DM
      integer, dimension(:,:), allocatable :: dmlist           !< List of neighbours for DM
      real(dblprec), dimension(:,:,:), allocatable :: dm_vect  !< Dzyaloshinskii-Moriya exchange vector
      ! Variables for Symmetric anisotropic exchange (SA)
      integer ::  max_no_saneigh                               !< Calculated maximum of neighbours for SA exchange
      integer, dimension(:), allocatable :: salistsize         !< Size of neighbour list for SA
      integer, dimension(:,:), allocatable :: salist           !< List of neighbours for SA
      real(dblprec), dimension(:,:,:), allocatable :: sa_vect  !< Symmetric anisotropic exchange vector
      ! Variables for CHIR exchange
      integer :: nn_chir_tot                                     !< Calculated number of neighbours with CHIR interactions
      integer ::  max_no_chirneigh                               !< Calculated maximum of neighbours for CHIR exchange
      integer, dimension(:), allocatable :: chirlistsize         !< Size of neighbour list for CHIR
      integer, dimension(:,:,:), allocatable :: chirlist           !< List of neighbours for CHIR
      real(dblprec), dimension(:,:), allocatable :: chir_coup  !< scalar chirality exchange coupling
      ! Variables for PD exchange
      integer :: nn_pd_tot                                     !< Calculated number of neighbours with PD interactions
      integer ::  max_no_pdneigh                               !< Calculated maximum of neighbours for PD exchange
      integer, dimension(:), allocatable :: pdlistsize         !< Size of neighbour list for PD
      integer, dimension(:,:), allocatable :: pdlist           !< List of neighbours for PD
      real(dblprec), dimension(:,:,:), allocatable :: pd_vect  !< Pseudo-Dipolar exchange vector
      ! Variables for BIQDM interactions
      integer :: nn_biqdm_tot                                     !< Calculated number of neighbours with BIQDM interactions
      integer ::  max_no_biqdmneigh                               !< Calculated maximum of neighbours for BIQDM exchange
      integer, dimension(:), allocatable :: biqdmlistsize         !< Size of neighbour list for BIQDM
      integer, dimension(:,:), allocatable :: biqdmlist           !< List of neighbours for BIQDM
      real(dblprec), dimension(:,:,:), allocatable :: biqdm_vect  !< BIQDM exchange vector
      ! Variables for BQ interactions
      integer :: nn_bq_tot                               !< Calculated number of neighbours with BQ interactions
      integer, dimension(:), allocatable :: bqlistsize   !< Size of neighbour list for BQ
      integer, dimension(:,:), allocatable :: bqlist     !< List of neighbours for BQ
      real(dblprec), dimension(:,:), allocatable :: j_bq !< Biquadratic exchange couplings
      ! Variables for four-spin ring (4SR) exchange interactions
      integer :: max_no_ringneigh                          !< Calculated maximum number of neighbours for 4SR exchange 
      integer :: nn_ring_tot                               !< Calculated number of neighbours with 4SR interactions
      integer, dimension(:), allocatable :: ringlistsize   !< Size of neighbour list for 4SR exchange
      integer, dimension(:,:,:), allocatable :: ringlist   !< List of neighbours for 4SR exchange
      real(dblprec), dimension(:,:), allocatable :: j_ring !< 4SR exchange couplings
      ! Variables for general four-spin exchange
      integer :: nn_fourx_tot                                     !< Calculated number of neighbours with CHIR interactions
      integer ::  max_no_fourxneigh                               !< Calculated maximum of neighbours for CHIR exchange
      integer, dimension(:), allocatable :: fourxlistsize         !< Size of neighbour list for CHIR
      integer, dimension(:,:,:), allocatable :: fourxlist           !< List of neighbours for CHIR
      real(dblprec), dimension(:,:), allocatable :: fourx_coup  !< scalar chirality exchange coupling
      ! Variables for anisotropy
      integer, dimension(:), allocatable :: taniso                !< Type of anisotropy (0-2)
      integer, dimension(:), allocatable :: taniso_diff           !< Type of anisotropy (0-2)
      real(dblprec), dimension(:), allocatable :: sb              !< Ratio between Cubic and Uniaxial anisotropy
      real(dblprec), dimension(:), allocatable :: sb_diff         !< Ratio between Cubic and Uniaxial anisotropy
      real(dblprec), dimension(:,:), allocatable :: kaniso        !< Anisotropy constant
      real(dblprec), dimension(:,:), allocatable :: kaniso_diff   !< Anisotropy constant
      real(dblprec), dimension(:,:), allocatable :: eaniso        !< Unit anisotropy vector
      real(dblprec), dimension(:,:), allocatable :: eaniso_diff   !< Unit anisotropy vector
      !!!!!!
      !!!!!!
      ! Variables for Biquadratic interaction 4spin-2site H11
      integer :: max_no_neigh_bqfull11									! Calculated maximum of neighbours for 4spin-2spin Biquadratic interactions (H11)
      integer, dimension(:), allocatable :: bqfull11listsize			! Size of neighbour list for 4spin-2site Biquadratic Interaction (H11)
      integer, dimension(:,:), allocatable :: bqfull11list				! Neighbours list for 4spin-2site Biquadratic Interaction (H11)
      real(dblprec), dimension(:,:,:,:), allocatable :: bqfull11_tens 	! Exchange tensor for 4spin-2site Biquadratic Interaction (H11)
      
      ! Variables for Biquadratic interaction 4spin-2site H21
      integer :: max_no_neigh_bqfull21									! Calculated maximum of neighbours for 4spin-2spin Biquadratic interactions (H21)
      integer, dimension(:), allocatable :: bqfull21listsize			! Size of neighbour list for 4spin-2site Biquadratic Interaction (H21)
      integer, dimension(:,:), allocatable :: bqfull21list				! Neighbours list for 4spin-2site Biquadratic Interaction (H21)
      real(dblprec), dimension(:,:), allocatable :: bqfull21			! Exchange scalar for 4spin-2site Biquadratic Interaction (H21)
      
      ! Variables for Biquadratic interaction 4spin-2site H22
      integer :: max_no_neigh_bqfull22									! Calculated maximum of neighbours for 4spin-2spin Biquadratic interactions (H22)
      integer, dimension(:), allocatable :: bqfull22listsize			! Size of neighbour list for 4spin-2site Biquadratic Interaction (H22)
      integer, dimension(:,:), allocatable :: bqfull22list				! Neighbours list for 4spin-2site Biquadratic Interaction (H22)
      real(dblprec), dimension(:,:,:,:), allocatable :: bqfull22_tens 	! Exchange tensor for 4spin-2site Biquadratic Interaction (H22)
      
      ! Variables for Biquadratic interaction 4spin-2site H23
      integer :: max_no_neigh_bqfull23									! Calculated maximum of neighbours for 4spin-2spin Biquadratic interactions (H23)
      integer, dimension(:), allocatable :: bqfull23listsize				! Size of neighbour list for 4spin-2site Biquadratic Interaction (H23)
      integer, dimension(:,:), allocatable :: bqfull23list				! Neighbours list for 4spin-2site Biquadratic Interaction (H23)
      real(dblprec), dimension(:,:,:), allocatable :: bqfull23_vec		! Exchange tensor for 4spin-2site Biquadratic Interaction (H23)

      ! Variables for Biquadratic interaction 4spin-2site H31
      integer :: max_no_neigh_bqfull31									! Calculated maximum of neighbours for 4spin-2spin Biquadratic interactions (H31)
      integer, dimension(:), allocatable :: bqfull31listsize			! Size of neighbour list for 4spin-2site Biquadratic Interaction (H31)
      integer, dimension(:,:), allocatable :: bqfull31list				! Neighbours list for 4spin-2site Biquadratic Interaction (H31)
      real(dblprec), dimension(:,:), allocatable :: bqfull31			! Exchange scalar for 4spin-2site Biquadratic Interaction (H31)

      ! Variables for Biquadratic interaction 4spin-2site H32
      integer :: max_no_neigh_bqfull32									! Calculated maximum of neighbours for 4spin-2spin Biquadratic interactions (H32)
      integer, dimension(:), allocatable :: bqfull32listsize			! Size of neighbour list for 4spin-2site Biquadratic Interaction (H32)
      integer, dimension(:,:), allocatable :: bqfull32list				! Neighbours list for 4spin-2site Biquadratic Interaction (H32)
      real(dblprec), dimension(:,:,:,:), allocatable :: bqfull32_tens 	! Exchange tensor for 4spin-2site Biquadratic Interaction (H32)

      ! Variables for Biquadratic interaction 4spin-2site H33
      integer :: max_no_neigh_bqfull33									! Calculated maximum of neighbours for 4spin-2spin Biquadratic interactions (H33)
      integer, dimension(:), allocatable :: bqfull33listsize			! Size of neighbour list for 4spin-2site Biquadratic Interaction (H33)
      integer, dimension(:,:), allocatable :: bqfull33list				! Neighbours list for 4spin-2site Biquadratic Interaction (H33)
      real(dblprec), dimension(:,:,:,:), allocatable :: bqfull33_tens 	! Exchange tensor for 4spin-2site Biquadratic Interaction (H33)

      ! Variables for Biquadratic interaction 4spin-2site H34
      integer :: max_no_neigh_bqfull34									! Calculated maximum of neighbours for 4spin-2spin Biquadratic interactions (H34)
      integer, dimension(:), allocatable :: bqfull34listsize			! Size of neighbour list for 4spin-2site Biquadratic Interaction (H34)
      integer, dimension(:,:), allocatable :: bqfull34list				! Neighbours list for 4spin-2site Biquadratic Interaction (H34)
      real(dblprec), dimension(:,:,:), allocatable :: bqfull34_vec		! Exchange tensor for 4spin-2site Biquadratic Interaction (H34)

      ! Variables for Biquadratic interaction 4spin-2site H35
      integer :: max_no_neigh_bqfull35										! Calculated maximum of neighbours for 4spin-2spin Biquadratic interactions (H35)
      integer, dimension(:), allocatable :: bqfull35listsize				! Size of neighbour list for 4spin-2site Biquadratic Interaction (H35)
      integer, dimension(:,:), allocatable :: bqfull35list					! Neighbours list for 4spin-2site Biquadratic Interaction (H35)
      real(dblprec), dimension(:,:,:,:,:), allocatable :: bqfull35_3tens	! Exchange tensor for 4spin-2site Biquadratic Interaction (H35)
      
      ! Variables for Biquadratic interaction 4spin-2site H36
      integer :: max_no_neigh_bqfull36										! Calculated maximum of neighbours for 4spin-2spin Biquadratic interactions (H36)
      integer, dimension(:), allocatable :: bqfull36listsize				! Size of neighbour list for 4spin-2site Biquadratic Interaction (H36)
      integer, dimension(:,:), allocatable :: bqfull36list					! Neighbours list for 4spin-2site Biquadratic Interaction (H36)
      real(dblprec), dimension(:,:,:,:,:), allocatable :: bqfull36_3tens	! Exchange tensor for 4spin-2site Biquadratic Interaction (H36)

      !!!!!!
      !!!!!!
      ! Variables for induced moments
      integer :: fix_num                                    !< Number of "fixed" moments
      integer :: max_no_neigh_ind                           !< Number of nearest neighbours for the induced treatment
      integer, dimension(:), allocatable :: fix_list        !< List containing the "fixed" moments
      integer, dimension(:), allocatable :: ind_list_full   !< Indication of whether a given moment is induced/fixed 1/0
      integer, dimension(:), allocatable :: ind_nlistsize   !< Size of the list for the induced moments
      integer, dimension(:,:), allocatable :: ind_nlist     !< Neighbour list between induced moments and their first permanent moments
      integer, dimension(:), allocatable :: fix_nlistsize   !< Size of the list for the permanent moments
      integer, dimension(:,:), allocatable :: fix_nlist     !< Neighbour list between permanent moments and their first induced moments
      real(dblprec), dimension(:), allocatable :: sus_ind   !< Scaling factor for the magneitc moment of the induced moments
      ! Variables for LSF
      integer, dimension(:), allocatable :: fs_nlistsize    !< Size of first shell neighbouring list for centered atom
      integer, dimension(:,:), allocatable :: nind          !< Index of firstshell-neighbour-list corresponds to neighbour-list
      integer, dimension(:,:), allocatable :: fs_nlist      !< First shell Neighbouring list for centered atom
      ! Variables for dipolar
      real(dblprec), dimension(:,:,:,:), allocatable :: Qdip         !< Matrix for dipole-dipole interaction
      real(dblprec), dimension(:,:,:,:), allocatable :: Qdip_macro   !< Matrix for macro spin dipole-dipole interaction
   end type ham_t

   public

end module HamiltonianDataType
