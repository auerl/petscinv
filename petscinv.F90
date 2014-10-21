!  ------------------------------------------------------------------------- 
!
!  PETScinv: 
!  
!  Solves a tomographic linear system in parallel with KSP.  Allows
!  for solution of rectangular systmes using a parallel implementaion
!  of lsqr, or the normal equations using various direct solvers
!
!  This is a complete Fortran 2008 rewrite of Lapo Boschi's original 
!  tomography toolbox. This is research code and we give no warranty
!  nor do we guaranty fitness for a particular purpose.
!
!  Compile with "make" and run with "./run_petscinv""
!
!  Copyright (c) 2013 - 2014 Ludwig Auer, ludwig.auer@tomo.ig.erdw.ethz.ch
!
!  ------------------------------------------------------------------------- 



!============================================================
! 
!  To specify all global parameters
!
      module global_param

        implicit none
        public


!***********************************
! petsc-headers for module_mesh

! petscsys.h - base PETSc routines      
! petscvec.h - vectors
! petscmat.h - matrices
! petscksp.h - Krylov subspace methods  
! petscpc.h  - preconditioners

#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscmat.h>
#include <finclude/petscis.h>
#include <finclude/petscao.h>
#include <finclude/petscpc.h>
#include <finclude/petscksp.h>
#include <finclude/petscviewer.h>
#include <finclude/petscdraw.h>
!***********************************

        real(8), parameter :: pi = 3.1415926535898D0
        real(8), parameter :: deg2rad = pi / 180.d0
        real(8), parameter :: rad2deg = 180.d0 / pi
        real(8), parameter :: r_earth = 6371.d0
        real(8), parameter :: deg2km = 111.32      ! 1 deg in km at the equator (approx.)
        real(8), parameter :: ref_eqinc = 5.0      ! reference grid, largest block size

        integer, parameter :: sp = selected_real_kind(6, 37)
        integer, parameter :: dp = selected_real_kind(15, 307)
        integer, parameter :: qp = selected_real_kind(33, 4931)        

        integer, parameter :: verbosity = 1        ! 1=standard, 2=a lot of output
        integer, parameter :: nlaym = 50           ! maximum number of layers
        integer, parameter :: n0max = 50000        ! maximum number of blocks per layer
        integer, parameter :: n1max = n0max*nlaym  ! maximum number of blocks per param
        integer, parameter :: npmax = 4            ! maximum number of physical params
        integer, parameter :: nlmax = 288          ! maximum number of latitude zones       
        integer, parameter :: matmx = 500          ! maximal number of matrices to read
        integer, parameter :: dmpmx = 500          ! maximal number of damping loops
        integer, parameter :: nprmx = 10000        ! maximal number per row for temporary arrays
        integer, parameter :: datmx = 5*(10**6)    ! maximal number of data per submatrix
        integer, parameter :: nperd = 30           ! maximal number of nonz in damping row

        PetscErrorCode   ierr
        PetscBool        flg
        PetscLogDouble   memory
        PetscInt         irun
        PetscInt         imat

        contains          

    !========================================================
    !
    !   converts an integer to a string
    !
          character(len=20) function int2str(k)
            integer, intent(in) :: k
            write(int2str, *) k
            int2str = adjustl(int2str)
          end function int2str
    !========================================================


      end module global_param
!============================================================



!============================================================
! 
!  module related to reading options from command line
!
      module module_options

        use global_param
        implicit none

        type,public :: opts
           PetscScalar    eqinc
           PetscScalar    eqinc_ref
           PetscInt       npars
           PetscInt       nlays
           PetscBool      adapt
           PetscBool      gvarr
           PetscChar(256) type
           PetscChar(256) sched   
           PetscChar(256) grdin                              
           PetscChar(256) inpar
           PetscChar(256) param
           PetscChar(256) synth
           PetscChar(256) refmod
           PetscChar(256) format
           PetscChar(256) projid
           character(5) :: parameter(4) ! maximum of 4 parameters
         contains
           procedure :: parse_options
           procedure :: tokenize_param
        end type opts

        contains

    !========================================================
    !
    !   reads options from command line
    !
          subroutine parse_options(this)

            implicit none
            class(opts) :: this

            ! Initialize some default variables
            this%eqinc_ref=ref_eqinc ! defined in the very top
            this%synth=''

            call PetscOptionsGetInt(PETSC_NULL_CHARACTER,&
                 '-number_of_layers',this%nlays,flg,ierr)
            call PetscOptionsGetReal(PETSC_NULL_CHARACTER,&
                 '-equatorial_increment',this%eqinc,flg,ierr)
            call PetscOptionsGetString(PETSC_NULL_CHARACTER,&
                 '-inversion_parameters',this%param,flg,ierr)
            call PetscOptionsGetString(PETSC_NULL_CHARACTER,&
                 '-reference_model',this%refmod,flg,ierr)
            call PetscOptionsHasName(PETSC_NULL_CHARACTER,&
                 '-adaptive_grid',this%adapt,ierr)
            call PetscOptionsHasName(PETSC_NULL_CHARACTER,&
                 '-grouped_varr',this%gvarr,ierr)
            call PetscOptionsGetString(PETSC_NULL_CHARACTER,&
                 '-matrix_schedule',this%sched,flg,ierr)
            call PetscOptionsGetString(PETSC_NULL_CHARACTER,&
                 '-output_format',this%format,flg,ierr)
            call PetscOptionsGetString(PETSC_NULL_CHARACTER,&
                 '-project_id',this%projid,flg,ierr)
            call PetscOptionsGetString(PETSC_NULL_CHARACTER,&
                 '-inparam_file',this%inpar,flg,ierr)
            call PetscOptionsGetString(PETSC_NULL_CHARACTER,&
                 '-synthetic_model',this%synth,flg,ierr)
            call PetscOptionsGetString(PETSC_NULL_CHARACTER,&
                 '-grid_info',this%grdin,flg,ierr)
            call PetscOptionsGetString(PETSC_NULL_CHARACTER,&
                 '-solution_type',this%type,flg,ierr)

          end subroutine parse_options


    !========================================================
    !
    !  identifies with which parameterization we are dealing with
    !
          subroutine tokenize_param(this)
            
            implicit none
            class(opts) :: this

            PetscInt pos1
            PetscInt pos2
            PetscInt i

            pos1 = 1
            this%npars = 0
 
            do
               pos2 = index(this%param(pos1:), ",")
               if (pos2 == 0) then
                  this%npars = this%npars + 1
                  this%parameter(this%npars) = this%param(pos1:)
                  exit
               end if
               this%npars = this%npars + 1
               this%parameter(this%npars) = this%param(pos1:pos1+pos2-2)
               pos1 = pos2+pos1
            end do
            
          end subroutine tokenize_param
    !========================================================

            


      end module module_options

!============================================================




!============================================================
! 
!  The module containing everything related to the mesh
!
      module module_mesh

        use global_param
        use module_options

        implicit none

        type,public :: mesh
           PetscInt       nlays
           PetscScalar    eqinc
           PetscBool      adapt
           PetscChar(256) grdin
           PetscScalar    eqinc_ref
           PetscInt       nlatzones

           ! To be defined
           PetscInt       blocks_per_param ! number of blocks per parameter
           PetscInt       blocks_all_param ! number of blocks per parameter

           PetscInt,      allocatable :: blocks_per_layer(:) ! number of blocks in each layer
           PetscInt,      allocatable :: nsqrs(:)
           PetscInt,      allocatable :: nsqrs_tot(:)
           
           PetscScalar,   allocatable :: xlamin(:,:)
           PetscScalar,   allocatable :: xlamax(:,:)
           PetscScalar,   allocatable :: xlomin(:,:)
           PetscScalar,   allocatable :: xlomax(:,:)

           PetscScalar,   allocatable :: locent(:,:)
           PetscScalar,   allocatable :: lacent(:,:)

           PetscScalar,   allocatable :: radmin(:,:)
           PetscScalar,   allocatable :: radmax(:,:)
           PetscInt,      allocatable :: levels(:,:)
           
         contains
           procedure :: setup_mesh
           procedure :: read_mesh
           procedure :: gen_mesh
           procedure :: coordinates => get_coordinates
        end type mesh

      contains

    !========================================================
    !
    !   procedure to initialize a new mesh
    !
        subroutine setup_mesh(this,options)
          
          implicit none
          class(mesh) :: this
          class(opts) :: options

          ! Local variables
          PetscBool                       exist
          PetscInt                        i,j          

          ! Allocate memory
          allocate ( this%blocks_per_layer( options%nlays ) )
          allocate ( this%levels( n0max , options%nlays ) )
          allocate ( this%xlamin( n0max , options%nlays ) )
          allocate ( this%xlamax( n0max , options%nlays ) )
          allocate ( this%xlomin( n0max , options%nlays ) )
          allocate ( this%xlomax( n0max , options%nlays ) )
          allocate ( this%radmin( n0max , options%nlays ) )
          allocate ( this%radmax( n0max , options%nlays ) )
          allocate ( this%locent( n0max , options%nlays ) )
          allocate ( this%lacent( n0max , options%nlays ) )
          allocate ( this%nsqrs_tot( nlmax ) ) ! 180/0.625
          allocate ( this%nsqrs( nlmax ) ) ! 180/0.625
          
          ! Initialize more mesh params
          this%eqinc = options%eqinc
          this%adapt = options%adapt
          this%nlays = options%nlays
          this%eqinc_ref = options%eqinc_ref
          this%blocks_per_layer(:) = 0
          this%blocks_per_param    = 0

          ! Select between variable and regular mesh
          select case ( options%adapt )
          case (.true.)

             call PetscPrintf(PETSC_COMM_WORLD,&
                  "MESHER: Reading mesh from disk\n",ierr)
             inquire ( file=options%grdin, exist=exist)
             if (.not. exist) then
                call PetscPrintf(PETSC_COMM_WORLD,&
                     "**** MESHER CRASHED **** : You didn't\
                      provide me with a grid info file!\n",ierr)
                stop
             end if
             call this%read_mesh(options%grdin)

          case (.false.)
             ! Generate a new grid based on eqincr
             call PetscPrintf(PETSC_COMM_WORLD,&
                  "MESHER: Generating mesh on the fly!\n",ierr)
             if ( options%eqinc.le.0.d5 ) then
                call PetscPrintf(PETSC_COMM_WORLD,&
                     "**** MESHER CRASHED ****: You didn't\
                      provide me with a valid eqincr value!\n",ierr)
                stop
             end if                          
             call this%gen_mesh()
          end select

          call PetscPrintf(PETSC_COMM_WORLD,"MESHER: Matrices are computed for "//&
               trim(int2str(options%npars))//" parameters\n",ierr)

          this%blocks_all_param = this%blocks_per_param*&
                                  options%npars   
                    
        end subroutine setup_mesh
    !========================================================



    !========================================================
    !
    !   subroutine to read a mesh from ascii on disk
    !
        subroutine read_mesh(this, grdin)

          implicit none
          class(mesh) :: this

          PetscChar(256)  grdin          
          PetscChar(5)    dummy
          PetscInt        i,j


          open(10,file=trim(grdin))           ! read levels of blocks
          open(20,file=trim(grdin)//".sh")    ! read coordinate range
          open(30,file=trim(grdin)//".lay")   ! read blocks per layer

          ! Read a stored grid from disk
          do i = 1,this%nlays
             read(30,*), this%blocks_per_layer(i)
             this%blocks_per_param = this%blocks_per_param + &
                                    this%blocks_per_layer(i)

             do j = 1,this%blocks_per_layer(i)
                read(10,*),dummy,this%levels(j,i)
                read(20,*),dummy,this%xlamin(j,i),&
                                    this%xlamax(j,i),&
                                    this%xlomin(j,i),&
                                    this%xlomax(j,i)

                ! Compute center lat and lon
                this%lacent(j,i)=abs((this%xlamax(j,i)-&
                     this%xlamin(j,i))/2) + this%xlamin(j,i)
                this%locent(j,i)=abs((this%xlomax(j,i)-&
                     this%xlomin(j,i))/2) + this%xlomin(j,i)                
             end do
          end do
             
          close(10)
          close(20)
          close(30)

          call PetscPrintf(PETSC_COMM_WORLD,&
               "MESHER: Mesh read successfully!\n",ierr)
          
        end subroutine read_mesh
    !========================================================



    !========================================================
    !
    !   subroutine to generate mesh on the fly
    !
        subroutine gen_mesh(this) !eqinc, eqinc_ref, nlays)
          
          implicit none
          class(mesh) :: this

          PetscScalar   coords( 6 )
          PetscInt      nsqrs_ref( nlmax )
          PetscScalar   colat, theta
          PetscScalar   delon, fact
          PetscInt      nlatz, nlatz_ref
          PetscInt      k, k_ref, n1lay
          PetscInt      i, j

          nlatz_ref = int( 180. / this%eqinc_ref )
          nlatz     = int( 180. / this%eqinc )
          colat     = -this%eqinc_ref / 2. 

          ! loop over latitudinal zones
          do k = 1,nlatz_ref             
             colat = colat + this%eqinc_ref
             theta = ( colat/180. )*pi
             delon = this%eqinc_ref / ( sin(theta) )
             nsqrs_ref(k) = int( 360./delon ) + 1
             if (mod(nsqrs_ref(k),2).ne.0) nsqrs_ref(k) = nsqrs_ref(k)-1
             ! correct nsqrs(k) such to be compatible with reference grid             
             if (360./nsqrs_ref(k).ge.this%eqinc_ref) then
                do while (mod(360./nsqrs_ref(k),this%eqinc_ref).ne.0)
                   nsqrs_ref(k) = nsqrs_ref(k) + 1
                end do
             elseif (360./nsqrs_ref(k).lt.this%eqinc_ref) then
                do while (mod(this%eqinc_ref,360./nsqrs_ref(k)).ne.0)
                   nsqrs_ref(k) = nsqrs_ref(k) - 1
                end do
             endif                              
             if (mod(nsqrs_ref(k),2).ne.0) then 
                call PetscPrintf(PETSC_COMM_WORLD,&
                "**** MESHER CRASHED ****: # pixels per\
                latitude has to be even\n",ierr)
                stop 
             end if
          enddo


          ! define actual grid as integer of reference grid
          fact  = this%eqinc_ref / this%eqinc
          n1lay = 0
          call PetscPrintf(PETSC_COMM_WORLD,&
               "MESHER: Ratio of coarse to fine grid: "//&
               trim(int2str(int(fact)))//"\n",ierr)
          do k = 1,nlatz
             k_ref = ( (k-1) / int(fact) ) + 1
             this%nsqrs(k)     = nsqrs_ref(k_ref) * int(fact)
             n1lay             = n1lay + this%nsqrs(k)
             this%nsqrs_tot(k) = n1lay
          enddo
          
          this%nlatzones = nlatz
          
          ! in all layers we have the same nr of blocks
          this%blocks_per_layer(:) = n1lay           
          
          call PetscPrintf(PETSC_COMM_WORLD,&
               'MESHER: Number of pixel with finest\
          parameterization: '//trim(int2str(&
               this%blocks_per_layer(1)))//'\n',ierr)
          
          ! continue down to other layers
          do i=1,this%nlays
             do j=1,this%blocks_per_layer(i)

                ! Get coordinate ranges
                coords = this%coordinates(j) 

                ! Minimum and maximum lat and lon
                this%levels(j,i) = 1 ! we just have one level
                this%xlamin(j,i) = coords(1)
                this%xlamax(j,i) = coords(2)
                this%xlomin(j,i) = coords(3)
                this%xlomax(j,i) = coords(4) ! dont need coords(5:6)

                ! Compute center lat and lon
                this%lacent(j,i) = coords(5)
                this%locent(j,i) = coords(6)
                
             end do
             this%blocks_per_param = this%blocks_per_param + &
                  this%blocks_per_layer(i)
          end do
          
        end subroutine gen_mesh
    !========================================================



    !========================================================
    !
    !    procedure to get center and range of a block
    !
        function get_coordinates(this,isqr) result(coords)

          implicit none
          class(mesh), intent(in) :: this
          PetscInt,    intent(in) :: isqr

          PetscScalar  coords(6)
          PetscScalar  xlamax,xlamin
          PetscScalar  xlomax,xlomin
          PetscScalar  xlomid,xlamid
          
          PetscScalar  rlong,rlati,rinlo
          PetscInt     ila,isq,ntot

          ntot   = 0
          xlomid = 0
          xlomax = 0
          xlomin = 0
          xlamid = 0
          xlamax = 0
          xlamin = 0

          ! loop(s) over all the blocks
          do ila=1,this%nlatzones
             ! increment latitude
             rlati=90.-(this%eqinc*(ila-1))
             ! calculate increment in longitude for this band
             rinlo=(360./this%nsqrs(ila))
             do isq=1,this%nsqrs(ila)
                rlong=(360./this%nsqrs(ila))*(isq-1)
                ntot=ntot+1
                if (ntot.eq.isqr) then
                   xlomid=rlong+(rinlo/2.d0)
                   xlamid=rlati-(this%eqinc/2.d0)
                   xlomax=rlong+(rinlo)
                   xlamax=rlati
                   xlomin=rlong
                   xlamin=rlati-(this%eqinc)
                   goto 600
                end if
             end do
          end do
          
600       continue ! gets here in case of top or bottom layer

          coords(1) = xlamin ! min lat
          coords(2) = xlamax ! max lat
          coords(3) = xlomin ! min lon
          coords(4) = xlomax ! max lon
          coords(5) = xlamid ! center lat
          coords(6) = xlomid ! center lon
          
        end function get_coordinates
    !========================================================


      end module module_mesh
!============================================================



!============================================================
! 
!  module related to the inversion schedule
!

      module module_schedule

        use module_mesh
        use module_options
        use global_param

        implicit none

        type,public :: sche
           PetscInt         mats_total ! number of submatrices 
           PetscInt         rows_total ! total number of rows
           PetscInt         loop_total
           PetscInt         nrdamprows
           PetscBool        irdamp
           PetscBool        indamp
           PetscBool        iddamp
           PetscBool        iscale
           PetscBool        isynth
           PetscChar(256)   inpar_file
           PetscChar(256)   sched_file
           PetscScalar,     allocatable :: ndamp(:,:)
           PetscScalar,     allocatable :: rdamp(:,:)
           PetscScalar,     allocatable :: ddamp(:,:)
           PetscScalar,     allocatable :: scale(:,:)
           PetscScalar,     allocatable :: layer(:)
           PetscScalar,     allocatable :: dloop(:,:)
           PetscInt,        allocatable :: prealloc(:) 
           PetscScalar,     allocatable :: meta_info(:,:)
           PetscInt,        allocatable :: fromto_info(:,:)
           character(500),  allocatable :: path_info(:,:)           
         contains
           procedure :: initialize_inversion
           procedure :: read_schedule_file
           procedure :: read_inparam_file
           procedure :: memory_prealloc
        end type sche


      contains

    !========================================================
    !
    !   procededure to initialize all inversion schedules
    !
        subroutine initialize_inversion(this,inopts,inmesh)

          implicit none

          class(sche) :: this
          class(opts) :: inopts
          class(mesh) :: inmesh

          this%inpar_file = inopts%inpar
          this%sched_file = inopts%sched

          this%loop_total = 0 ! number of global damping weights
          this%rows_total = 0 ! total number of measurements
          this%mats_total = 0 ! total number of submatrices
          this%nrdamprows = 0 ! default is 0

          ! logicals on the type of damping applied
          this%irdamp = .false.
          this%indamp = .false.
          this%iddamp = .false.
          this%iscale = .false.
          
          ! logicals about synthetic testing
          if (trim(inopts%synth)=='') then
             this%isynth = .false.          
          else
             this%isynth = .true.
          end if

          ! allocate submatrix and damping schedule
          allocate ( this%rdamp( inopts%nlays, npmax+1 ) )
          allocate ( this%ndamp( inopts%nlays, npmax ) )
          allocate ( this%ddamp( inopts%nlays, 2 ) )
          allocate ( this%scale( inopts%nlays, 3 ) )
          allocate ( this%layer( inopts%nlays+1 ) )

          allocate ( this%path_info( matmx, 4 ) )
          allocate ( this%meta_info( matmx, 4 ) )
          allocate ( this%fromto_info( matmx, 2 ) )

          allocate ( this%dloop( 4, dmpmx ) )

          ! read damping and matrix schedules
          call this % read_inparam_file ( inopts, inmesh )
          call this % read_schedule_file ( inopts )

          ! now we now how large the matrix is going to be
          ! so allocate space for preallocvec
          allocate ( this%prealloc ( this%nrdamprows + &
                                     this%rows_total ) )
          
          call this % memory_prealloc ( )
                   
        end subroutine initialize_inversion
    !========================================================



    !========================================================
    !
    !   Procededure to read a damping parameters and other
    !   necessary information from the inparam file
    !
        subroutine read_inparam_file(this,inopts,inmesh)

          implicit none

          class(sche):: this
          class(opts):: inopts
          class(mesh):: inmesh

          PetscInt  ios,j,i
          PetscInt, parameter :: fh=16 ! file handler
          PetscChar(256) line,key,val,dummy

          open(unit=fh,file=this%inpar_file,iostat=ios)
          ios = 0
          do
             read(fh, fmt='(a256)', iostat=ios) line               

             if (ios < 0) then
                call PetscPrintf(PETSC_COMM_WORLD,&
                     "SCHEDULER: Read damping scheme successfully\n",ierr)
                exit
             end if

             if (len(trim(line)) < 1 .or.&
                   line(1:1) == ';') then
                continue
             else              
                 select case(trim(line))              
                   case ('RDAMP_DEPTH')
                      call PetscPrintf(PETSC_COMM_WORLD,&
                           "SCHEDULER: Reading roughness damping scheme\n",ierr)
                      read(fh, fmt=*, iostat=ios) dummy
                      do i=1,inopts%nlays
                         read(fh, fmt=*, iostat=ios) this%rdamp(i,1:5)
                         if ( verbosity > 1 ) print*,this%rdamp(i,:)
                      end do
                   case ('NDAMP_DEPTH')
                      call PetscPrintf(PETSC_COMM_WORLD,&
                           "SCHEDULER: Reading norm damping scheme\n",ierr)
                      read(fh, fmt=*, iostat=ios) dummy
                      do i=1,inopts%nlays
                         read(fh, fmt=*, iostat=ios) this%ndamp(i,1:4)
                         if ( verbosity > 1 ) print*,this%ndamp(i,:)
                      end do
                   case ('DDAMP_DEPTH')
                      call PetscPrintf(PETSC_COMM_WORLD,&
                           "SCHEDULER: Reading difference damping scheme\n",ierr)
                      read(fh, fmt=*, iostat=ios) dummy
                      do i=1,inopts%nlays
                         read(fh, fmt=*, iostat=ios) this%ddamp(i,1:2)
                         if ( verbosity > 1 ) print*,this%ddamp(i,1:2)
                      end do
                   case ('SCALING_COEFFS')
                      call PetscPrintf(PETSC_COMM_WORLD,&
                           "SCHEDULER: Reading vp-to-vs scaling scheme\n",ierr)
                      read(fh, fmt=*, iostat=ios) dummy
                      do i=1,inopts%nlays
                         read(fh, fmt=*, iostat=ios) this%scale(i,1:3)                         
                         if ( verbosity > 1 ) print*,this%scale(i,1:3)
                      end do
                   case ('LAYER_BOTTOMS')
                      call PetscPrintf(PETSC_COMM_WORLD,&
                           "SCHEDULER: Reading layer bottom depths for postprocessing\n",ierr)
                      read(fh, fmt=*, iostat=ios) dummy
                      this%layer(1)=0.d0 ! 0th index interface always surface = depth zero
                      do i=2,inopts%nlays+1
                         read(fh, fmt=*, iostat=ios) this%layer(i)
                         if ( verbosity > 1 ) print*,this%layer(i)
                      end do
                   case ('DAMP_LOOP')
                      call PetscPrintf(PETSC_COMM_WORLD,&
                           "SCHEDULER: Reading global damping loop scheme\n",ierr)
                      read(fh, fmt=*, iostat=ios) dummy
                      read(fh, fmt=*, iostat=ios) this%loop_total
                      do i=1,4
                         read(fh, fmt=*, iostat=ios) &
                              this%dloop(i,1:this%loop_total)
                         if ( verbosity > 1 ) print*,this%dloop(i,1:this%loop_total)
                      end do
                      ! determine total number of damping rows
                      do i=1,this%loop_total
                         ! nr of roughness damping rows
                         if ( this%dloop(1,i) > 0. ) this%irdamp = .true.
                         if ( this%dloop(2,i) > 0. ) this%indamp = .true.
                         if ( this%dloop(3,i) > 0. ) this%iddamp = .true.
                         if ( this%dloop(4,i) > 0. ) this%iscale = .true.
                      end do
                      if ( this%indamp ) then
                         this%nrdamprows = this%nrdamprows + &
                              inmesh%blocks_per_param * &
                              inopts%npars
                      end if
                      if ( this%irdamp ) then
                         this%nrdamprows = this%nrdamprows + &
                              inmesh%blocks_per_param * &
                              inopts%npars
                      end if
                      if ( this%iddamp ) then
                         if (inopts%npars==4) then
                            this%nrdamprows = this%nrdamprows + &
                                 inmesh%blocks_per_param*2
                         else if ((inopts%npars==2).or.(inopts%npars==3)) then
                            this%nrdamprows = this%nrdamprows + &
                                 inmesh%blocks_per_param
                         end if
                      end if                      
                      if ( this%iscale ) then
                         if (inopts%npars==4) then
                            this%nrdamprows = this%nrdamprows + &
                                 inmesh%blocks_per_param*2
                         else if (inopts%npars==2) then
                            this%nrdamprows = this%nrdamprows + &
                                 inmesh%blocks_per_param
                         end if
                      end if                      
                      call PetscPrintf(PETSC_COMM_WORLD,"SCHEDULER: In total we attach : "//&
                           trim(int2str(this%nrdamprows))//" damping rows\n",ierr)
                   end select
                end if
             end do
             close(fh)

        end subroutine read_inparam_file
    !========================================================



    !========================================================
    !
    !   procededure to read a matrix schedule
    !
        subroutine read_schedule_file(this,inopts)

          implicit none

          class(sche):: this          
          class(opts):: inopts

          PetscInt  ios,j,i
          PetscInt, parameter :: fh=15 ! file handler

          ios = 0
          open(fh,file=trim(this%sched_file))
          this%mats_total = 1

          do
             call PetscPrintf(PETSC_COMM_WORLD,&
                  "SCHEDULER: Reading meta info for submatrix "//&
                   trim(int2str(this%mats_total))//"\n",ierr)
             do j=1,4
                read(unit=fh,fmt=*,iostat=ios) &
                     this%path_info( this%mats_total, j )
                 if (verbosity >1 ) print*,trim(this%path_info( this%mats_total, j ))
             end do
             do j=1,4
                read(unit=fh,fmt=*,iostat=ios) &
                     this%meta_info( this%mats_total, j )
                if (verbosity > 1 ) print*,this%meta_info( this%mats_total, j )
             end do
             if ( trim(this%path_info( this%mats_total, 1 )) =='EOF' ) then
                ! remove trash written in the last entry
                this%meta_info( this%mats_total, : ) = 0.
                this%path_info( this%mats_total, : ) = ''
                this%mats_total = this%mats_total - 1
                exit
             else                
                call PetscPrintf(PETSC_COMM_WORLD,"SCHEDULER: Subset contains  "//&
                trim(int2str(int(this%meta_info(this%mats_total,4))))//" rows\n",ierr)
                this%rows_total = this%rows_total + &
                     int(this%meta_info( this%mats_total, 4 ))
                this%mats_total = this%mats_total + 1
             end if
          end do
         close(fh)
                 
         call PetscPrintf(PETSC_COMM_WORLD,&
              "SCHEDULER: We will combine a total of: "//&
              trim(int2str(this%rows_total))//" data\n",ierr)
         call PetscPrintf(PETSC_COMM_WORLD,&
              "SCHEDULER: Successfully read submatrix schedule!\n ",ierr)

       end subroutine read_schedule_file
    !========================================================       



    !========================================================
    !
    !   procededure to define
    !
        subroutine memory_prealloc(this)

          implicit none

          class(sche):: this
         
          PetscInt i,j,row_glob
          PetscInt, allocatable :: pnt(:)
                  
          ! loop over all submatrices, detailed preallocation
          row_glob = 0
          do i=1,this%mats_total
             ! meta_info(i,3) contains nr rows of submatrix
             allocate( pnt ( int(this%meta_info(i,4)) ) )
             ! path_info(i,3) contains pointer vector
             open(70,file=trim(this%path_info(i,3)),status='old')     
             ! loop over rows in submatrix
             do j=1,int(this%meta_info(i,4))
                row_glob=row_glob+1
                read(70,*) pnt(j)
                this%prealloc(row_glob) = pnt(j)-pnt(j-1) 
             end do
             deallocate ( pnt )
             close(70)
          end do
          ! and now the same for all damping rows
          do i=1,this%nrdamprows
             row_glob=row_glob+1
             ! max nr of entries per damping row
             this%prealloc(row_glob) = nperd 
          end do

         call PetscPrintf(PETSC_COMM_WORLD,"SCHEDULER: Successfully"//&
         "pre-allocated memory!\n ",ierr)

        end subroutine memory_prealloc
    !========================================================       


      end module module_schedule
!============================================================



!============================================================
! 
!  module containing everything related to reading matrices
!
      module module_matrix

          use global_param
      
          use module_schedule
          use module_options
          use module_mesh
 
          implicit none

          type,public :: matr
             Mat          A   ! global system matrix
             Vec          x   ! solution vector
             Vec          b   ! rhs vector
             Vec          b_store   ! to store the rhs vector
             Vec          b_synth   ! synthetics rhs vector
             Vec          b_adotx   ! rhs vector for cummulative variance red
             Vec          b_dummy   ! rhs vector for grouped variance red
             Vec          b_rough   ! synthetics rhs vector
             Vec          x_synth   ! synthetics model vector
             Vec          mask_dmp   ! masking the damping rows
             Vec          mask_sub   ! masking everything but subset rows
             Vec          mask_dat   ! masking everything but roughness damping rows
             Mat          ATb ! A transposed * b
             Mat          ATA ! bormal equations
             PetscInt     row
             PetscInt     row_start
             PetscInt     row_end      
             PetscMPIInt  processor
           contains
             procedure :: initialize_matrix
             procedure :: read_submatrices
             procedure :: assemble_matrix
             procedure :: assemble_vectors
             procedure :: apply_rdamp
             procedure :: apply_ddamp
             procedure :: apply_ndamp
             procedure :: apply_scale
             procedure :: store_rhs
             procedure :: restore_rhs
             procedure :: compute_adotx
             procedure :: destroy_matrix
             procedure :: read_synth_model
          end type matr

        contains

    !========================================================
    !
    !   procededure to initialize linear system to solve
    !
          subroutine initialize_matrix(this,insche,inmesh)

            implicit none

            class(matr) :: this
            class(sche) :: insche
            class(mesh) :: inmesh

            ! Get the rank of this processor
            call MPI_Comm_rank ( PETSC_COMM_WORLD, this%processor, ierr)

            call MatCreate ( PETSC_COMM_WORLD, this%A, ierr )

            ! changing size of A matrix, to make space for damping rows
            call MatSetSizes   ( this%A, PETSC_DECIDE,             &
                                 PETSC_DECIDE, insche%rows_total+ &
                                 insche%nrdamprows, &
                                 inmesh%blocks_all_param, ierr)

            ! set type of damping matrix
            call MatSetType    ( this%A, MATAIJ, ierr)

            ! Need to set up first before ownership
            call MatSetUp ( this%A, ierr)
            
            ! Get owenership range for this processor
            call MatGetOwnershipRange( this%A, this%row_start, &
                                       this%row_end, ierr)

            ! Preallocate space for more efficient matrix assembly
            call MatMPIAIJSetPreallocation ( this%A, PETSC_NULL_INTEGER, & 
                 insche%prealloc(this%row_start+1:this%row_end), &
                 PETSC_NULL_INTEGER, &
                 insche%prealloc(this%row_start+1:this%row_end), ierr )
            call MatSeqAIJSetPreallocation ( this%A, PETSC_NULL_INTEGER, &
                 insche%prealloc(this%row_start+1:this%row_end), ierr)
            
            if (verbosity > 1) print*,"Processor:",this%processor
            if (verbosity > 1) print*,this%row_start,this%row_end

            ! Create the rhs vectors
            call VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, &
                                insche%rows_total+ &
                                insche%nrdamprows, &
                                this%b,ierr)
            call VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, &
                                insche%rows_total+ &
                                insche%nrdamprows, &
                                this%b_synth,ierr)
            call VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, &
                                insche%rows_total+ &
                                insche%nrdamprows, &
                                this%b_store,ierr)
            call VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, &
                                insche%rows_total+ &
                                insche%nrdamprows, &
                                this%b_adotx,ierr)
            call VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, &
                                insche%rows_total+ &
                                insche%nrdamprows, &
                                this%b_dummy,ierr)
            call VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, &
                                insche%rows_total+ &
                                insche%nrdamprows, &
                                this%b_rough,ierr)

            ! Mask vectors for VARR, ROUGH and SYNTH
            call VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, &
                                insche%rows_total+ &
                                insche%nrdamprows, &
                                this%mask_dmp,ierr)
            call VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, &
                                insche%rows_total+ &
                                insche%nrdamprows, &
                                this%mask_dat,ierr)
            call VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, &
                                insche%rows_total+ &
                                insche%nrdamprows, &
                                this%mask_sub,ierr)

            ! Solution vectors
            call VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, &
                                inmesh%blocks_all_param, this%x, ierr)
            call VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, &
                                inmesh%blocks_all_param, this%x_synth, ierr)

            call PetscPrintf(PETSC_COMM_WORLD,"MATRIX: Initialized global "//&
            "kernel matrix and the rhs vectors\n",ierr)

          end subroutine initialize_matrix
    !========================================================



    !========================================================
    !
    !   procedure to determine the damping operator
    !
          subroutine apply_rdamp(this,inopts,inmesh,insche,irun) 
            
            implicit none

            class(matr) :: this
            class(opts) :: inopts
            class(mesh) :: inmesh
            class(sche) :: insche

            PetscScalar,      dimension(4) :: weight
            PetscInt,   dimension(6,nperd) :: ib
            PetscInt,         dimension(6) :: kk

            PetscInt      p,z,i,j,k,u,v
            PetscInt      nc,zo,ia,irun
            PetscScalar   weightv

            PetscInt      numlay,nparm,clay
            PetscScalar   lomin,lomax,lamin,lamax
            PetscScalar   romin,romax,ramin,ramax

            PetscInt      ind_petsc
            PetscScalar   val_petsc
            PetscScalar   rhs_petsc

            if (verbosity > 1) call PetscPrintf(PETSC_COMM_WORLD,&
                 "Starting with RDAMP, at row "//&
                 trim(int2str(this%row))//"\n",ierr)

            if ( this%processor == 0 ) then ! only estimate rdamp operator with proc 0          
               
               ! Loop over inversion parameters
               do u=1,inopts%npars
                  if (verbosity > 1) print*,"rdamp, parameter: ",u
                  
                  ! Loop over layers
                  clay=0
                  do i=1,inopts%nlays
                     if (verbosity > 1) print*,"rdamp, layer: ",i
                     

                     if ( irun == 0 ) then
                        ! The special case in which we initialize 
                        ! the roughness operator with uniform weight
                        ! for L-curve roughness estimation
                        weight    = 1.d0
                        weightv   = 1.d0
                     else 
                        ! dampig weights at different blocks sizes
                        weight(1) = insche%rdamp(i,u) * 1.0 * insche%dloop(1,irun)
                        weight(2) = insche%rdamp(i,u) * 0.9 * insche%dloop(1,irun)
                        weight(3) = insche%rdamp(i,u) * 0.8 * insche%dloop(1,irun)
                        weight(4) = insche%rdamp(i,u) * 0.7 * insche%dloop(1,irun)
                        ! column 5 in damping scheme is vertical
                        weightv   = insche%rdamp(i,5)                        
                     end if
                                 
                     do j=1,inmesh%blocks_per_layer(i)
                        
                        ! identify the current pixel index
                        ia=(u-1)*inmesh%blocks_per_param+clay+j               

                        if (ia.gt.(inmesh%blocks_per_param*inopts%npars)) then                 
                           print*,"RDAMP CRASHED: Undefined block index!"
                           stop 
                        end if
                        
                        ! this speeds things up a bit
                        ramin=inmesh%xlamin(j,i)
                        ramax=inmesh%xlamax(j,i)
                        romin=inmesh%xlomin(j,i)
                        romax=inmesh%xlomax(j,i)
 
                        kk(:)=0
                        ib(:,:)=0
                                                                        
                        ! 2nd loop over all pixels in matrix
                        do z=1,3
                           if(z.eq.1)then
                              zo=i ! this is the curr
                           elseif(z.eq.2)then
                              if(i.eq.1)then
                                 goto 999 ! skip this top layer, we are in first layer
                              else
                                 zo=i-1
                              endif
                           elseif(z.eq.3)then
                              if(i.eq.inopts%nlays)then
                                 goto 999 ! skip this top layer, we are in first layer
                              else
                                 zo=i+1
                              endif
                           endif
                           
                           do k=1,inmesh%blocks_per_layer(zo)
                              
                              ! this speeds things up a bit
                              lamin = inmesh%xlamin(k,zo)
                              lamax = inmesh%xlamax(k,zo)
                              lomin = inmesh%xlomin(k,zo)
                              lomax = inmesh%xlomax(k,zo)
                              
                              if(z.eq.3)then
                                 if((lamin.ge.ramin).and.(lamax.le.ramax).and.&
                                   (lomin.ge.romin).and.(lomax.le.romax))then
                                    kk(5)=kk(5)+1
                                    ib(5,kk(5))=(u-1)*inmesh%blocks_per_param+&
                                         (clay+inmesh%blocks_per_layer(i))+k
                                 elseif((ramin.ge.lamin).and.(ramax.le.lamax).and.&
                                   (romin.ge.lomin).and.(romax.le.lomax))then
                                    kk(5)=kk(5)+1
                                    ib(5,kk(5))=(u-1)*inmesh%blocks_per_param+&
                                         (clay+inmesh%blocks_per_layer(i))+k
                                 endif
                              elseif(z.eq.2)then
                                 if((lamin.ge.ramin).and.(lamax.le.ramax).and.&
                                   (lomin.ge.romin).and.(lomax.le.romax))then
                                    kk(6)=kk(6)+1
                                    ib(6,kk(6))=(u-1)*inmesh%blocks_per_param+&
                                         (clay-inmesh%blocks_per_layer(i-1))+k
                                 elseif((ramin.ge.lamin).and.(ramax.le.lamax).and.&
                                   (romin.ge.lomin).and.(romax.le.lomax))then
                                    kk(6)=kk(6)+1
                                    ib(6,kk(6))=(u-1)*inmesh%blocks_per_param+&
                                         (clay-inmesh%blocks_per_layer(i-1))+k
                                 endif
                              elseif(z.eq.1)then
                                 if((lamin.eq.ramin).and.(lamax.eq.ramax).and.&
                                   (lomin.eq.romin).and.(lomax.eq.romax))then
                                    continue ! do nothing actually
                                 elseif(lamin.eq.ramax)then ! south adjacent
                                    if((lomin.ge.romin).and.(lomax.le.romax))then
                                       kk(1)=kk(1)+1
                                       ib(1,kk(1))=(u-1)*inmesh%blocks_per_param+clay+k
                                    elseif((romin.ge.lomin).and.(romax.le.lomax))then
                                       kk(1)=kk(1)+1
                                       ib(1,kk(1))=(u-1)*inmesh%blocks_per_param+clay+k
                                    elseif((lomin.lt.romax).and.(lomax.gt.romax))then
                                       kk(1)=kk(1)+1
                                       ib(1,kk(1))=(u-1)*inmesh%blocks_per_param+clay+k
                                    elseif((lomin.lt.romin).and.(lomax.gt.romin))then
                                       kk(1)=kk(1)+1
                                       ib(1,kk(1))=(u-1)*inmesh%blocks_per_param+clay+k
                                    endif
                                 elseif(lamax.eq.ramin)then ! north adjacent
                                    if((lomin.ge.romin).and.(lomax.le.romax))then
                                       kk(2)=kk(2)+1
                                       ib(2,kk(2))=(u-1)*inmesh%blocks_per_param+clay+k
                                    elseif((romin.ge.lomin).and.(romax.le.lomax))then
                                       kk(2)=kk(2)+1
                                       ib(2,kk(2))=(u-1)*inmesh%blocks_per_param+clay+k
                                    elseif((lomin.lt.romax).and.(lomax.gt.romax))then
                                       kk(2)=kk(2)+1
                                       ib(2,kk(2))=(u-1)*inmesh%blocks_per_param+clay+k
                                    elseif((lomin.lt.romin).and.(lomax.gt.romin))then
                                       kk(2)=kk(2)+1
                                       ib(2,kk(2))=(u-1)*inmesh%blocks_per_param+clay+k
                                    endif
                                 elseif(lomin.eq.romax)then ! east adjacent
                                    if((lamin.ge.ramin).and.(lamax.le.ramax))then
                                       kk(3)=kk(3)+1
                                       ib(3,kk(3))=(u-1)*inmesh%blocks_per_param+clay+k
                                    elseif((ramin.ge.lamin).and.(ramax.le.lamax))then
                                       kk(3)=kk(3)+1
                                       ib(3,kk(3))=(u-1)*inmesh%blocks_per_param+clay+k
                                    endif
                                 elseif(lomax.eq.romin)then ! west adjacent
                                    if((lamin.ge.ramin).and.(lamax.le.ramax))then
                                       kk(4)=kk(4)+1
                                       ib(4,kk(4))=(u-1)*inmesh%blocks_per_param+clay+k
                                    elseif((ramin.ge.lamin).and.(ramax.le.lamax))then
                                       kk(4)=kk(4)+1
                                       ib(4,kk(4))=(u-1)*inmesh%blocks_per_param+clay+k
                                    endif
                                    ! special treatment at boundary
                                 elseif((lomax.eq.360.).and.(romin.eq.0.))then ! west
                                    if((lamin.ge.ramin).and.(lamax.le.ramax))then
                                       kk(4)=kk(4)+1
                                       ib(4,kk(4))=(u-1)*inmesh%blocks_per_param+clay+k
                                    elseif((ramin.ge.lamin).and.(ramax.le.lamax))then
                                       kk(4)=kk(4)+1
                                       ib(4,kk(4))=(u-1)*inmesh%blocks_per_param+clay+k
                                    endif
                                 elseif((lomin.eq.0.).and.(romax.eq.360.))then ! east
                                    if((lamin.ge.ramin).and.(lamax.le.ramax))then
                                       kk(3)=kk(3)+1
                                       ib(3,kk(3))=(u-1)*inmesh%blocks_per_param+clay+k
                                    elseif((ramin.ge.lamin).and.(ramax.le.lamax))then
                                       kk(3)=kk(3)+1
                                       ib(3,kk(3))=(u-1)*inmesh%blocks_per_param+clay+k
                                    endif
                                 endif
                              endif
                           end do ! end of loop over pixels in same layer
999                        continue ! gets here in case of top or bottom layer
                        end do ! loop over top mid and bot layer

               
                        ind_petsc=ia-1 ! because petsc uses 0 based indices
                        val_petsc=0
                        
                        do p=1,kk(1)
                           val_petsc=val_petsc+real((1./real(kk(1)))*&
                                weight(inmesh%levels(j,i)))            
                        end do
                        do p=1,kk(2)
                           val_petsc=val_petsc+real((1./real(kk(2)))*&
                                weight(inmesh%levels(j,i)))                        
                        end do
                        do p=1,kk(3)
                           val_petsc=val_petsc+real((1./real(kk(3)))*&
                                weight(inmesh%levels(j,i)))            
                        end do
                        do p=1,kk(4)
                           val_petsc=val_petsc+real((1./real(kk(4)))*&
                                weight(inmesh%levels(j,i)))            
                        end do
                        do p=1,kk(5)
                           val_petsc=val_petsc+real((1./real(kk(5)))*&
                                weight(inmesh%levels(j,i))*weightv)            
                        end do
                        do p=1,kk(6)
                           val_petsc=val_petsc+real((1./real(kk(6)))*&
                                weight(inmesh%levels(j,i))*weightv)            
                        end do
                        
                        call MatSetValue(this%A,this%row,ind_petsc,val_petsc,INSERT_VALUES,ierr)
                  
                        ! Increment the matrix, off-diagonal entries
                        do v=1,kk(1)
                           val_petsc=real((-1./real(kk(1)))*weight(inmesh%levels(j,i)))
                           ind_petsc=ib(1,v)-1
                           call MatSetValue(this%A,this%row,ind_petsc,val_petsc,INSERT_VALUES,ierr)
                           if(ib(1,v).gt.(inmesh%blocks_per_param*inopts%npars))then
                              print*,"RDAMP CRASHED: undefined index"
                              stop
                           endif
                        end do
                        do v=1,kk(2)
                           val_petsc=real((-1./real(kk(2)))*weight(inmesh%levels(j,i)))
                           ind_petsc=ib(2,v)-1
                           call MatSetValue(this%A,this%row,ind_petsc,val_petsc,INSERT_VALUES,ierr)
                           if(ib(2,v).gt.(inmesh%blocks_per_param*inopts%npars))then
                              print*,"RDAMP CRASHED: undefined index"
                              stop
                           endif
                        end do
                        do v=1,kk(3)
                           val_petsc=real((-1./real(kk(3)))*weight(inmesh%levels(j,i)))
                           ind_petsc=ib(3,v)-1
                           call MatSetValue(this%A,this%row,ind_petsc,val_petsc,INSERT_VALUES,ierr)
                           if(ib(3,v).gt.(inmesh%blocks_per_param*inopts%npars))then
                              print*,"RDAMP CRASHED: undefined index"
                              stop
                           endif
                        end do
                        do v=1,kk(4)
                           val_petsc=real((-1./real(kk(4)))*weight(inmesh%levels(j,i)))
                           ind_petsc=ib(4,v)-1
                           call MatSetValue(this%A,this%row,ind_petsc,val_petsc,INSERT_VALUES,ierr)
                           if(ib(4,v).gt.(inmesh%blocks_per_param*inopts%npars))then
                              print*,"RDAMP CRASHED: undefined index"
                              stop
                           endif
                        end do
                        do v=1,kk(5)
                           val_petsc=real((-1./real(kk(5)))*weight(inmesh%levels(j,i))*weightv)
                           ind_petsc=ib(5,v)-1
                           call MatSetValue(this%A,this%row,ind_petsc,val_petsc,INSERT_VALUES,ierr)
                           if(ib(5,v).gt.(inmesh%blocks_per_param*inopts%npars))then
                              print*,"RDAMP CRASHED: undefined index"
                              stop
                           endif
                        end do
                        do v=1,kk(6)
                           val_petsc=real((-1./real(kk(6)))*weight(inmesh%levels(j,i))*weightv)
                           ind_petsc=ib(6,v)-1
                           call MatSetValue(this%A,this%row,ind_petsc,val_petsc,INSERT_VALUES,ierr)
                           if(ib(6,v).gt.(inmesh%blocks_per_param*inopts%npars))then
                              print*,"RDAMP CRASHED: undefined index"
                              stop
                           endif
                        end do
                        
                        ! increment rhs and pnt vectors
                        rhs_petsc=0.
                        call VecSetValue(this%b,this%row,rhs_petsc,INSERT_VALUES,ierr)
                        call VecSetValue(this%mask_dmp,this%row,rhs_petsc,INSERT_VALUES,ierr)
                        call VecSetValue(this%mask_dat,this%row,1.d0,INSERT_VALUES,ierr)
                        
                        this%row=this%row+1

                     end do
                     clay=clay+inmesh%blocks_per_layer(i) ! we go to next layer
                  end do
               end do
               
            end if ! end if rank
        
          end subroutine apply_rdamp
    !========================================================



    !========================================================
    !
    !   subroutine to append norm damping tems
    !
    !   NOTE: This routine has just been adapted from an old 
    !         code by Lapo and needs to be checked & tested
    !         
    !
          subroutine apply_ndamp(this,inopts,inmesh,insche,irun) 
            
            implicit none

            class(matr) :: this
            class(opts) :: inopts
            class(mesh) :: inmesh
            class(sche) :: insche

            PetscScalar   val,rhs
            PetscInt      ind,irun
            PetscInt      u,i,j
           
            if (verbosity > 1) call PetscPrintf(PETSC_COMM_WORLD,&
                 "Starting with NDAMP, at row "//&
                 trim(int2str(this%row))//"\n",ierr)

            if ( this%processor == 0 ) then ! only estimate rdamp operator with proc 0              
               ind=0
               rhs=0. ! always 0
               do u=1,inopts%npars
                  if (verbosity > 1) print*,"ndamp, parameter: ",u
                  do i=1,inopts%nlays
                     if (verbosity > 1) print*,"ndamp, layer: ",i                  
                     do j=1,inmesh%blocks_per_layer(i)
                        ! contains a scaling with the area of a voxel
                        val=insche%ndamp(i,u)*insche%dloop(2,irun)*&
                             (inmesh%eqinc_ref/real(inmesh%levels(j,i)))**2
                        call MatSetValue(this%A,this%row,ind,val,INSERT_VALUES,ierr)
                        call VecSetValue(this%b,this%row,rhs,INSERT_VALUES,ierr)
                        call VecSetValue(this%mask_dmp,this%row,rhs,INSERT_VALUES,ierr)
                        call VecSetValue(this%mask_dat,this%row,rhs,INSERT_VALUES,ierr)
                        this%row=this%row+1 
                        ind=ind+1
                     end do
                  end do
               end do                                                                    
            end if

          end subroutine apply_ndamp
    !========================================================



    !========================================================
    !
    !   procedure to determine the damping operator
    !
          subroutine apply_ddamp(this,inopts,inmesh,insche,irun) 
            
            implicit none

            class(matr) :: this
            class(opts) :: inopts
            class(mesh) :: inmesh
            class(sche) :: insche

            PetscScalar   val,rhs
            PetscInt      ind,irun
            PetscInt      i,j,k

            k=0
            rhs=0.

            if (verbosity > 1) call PetscPrintf(PETSC_COMM_WORLD,&
                 "Starting with DDAMP, at row "//&
                 trim(int2str(this%row))//"\n",ierr)

            if ( this%processor == 0 ) then ! only estimate rdamp operator with proc 0              
               do i=1,inopts%nlays

                  if (verbosity > 1) print*,"ddamp, layer: ",i
                  do j=k+1,k+inmesh%blocks_per_layer(i) ! normal indices
                     
                     ! damp the difference between parameters 2 and 1
                     val=insche%ddamp(i,1)*insche%dloop(3,irun)
                     ind=j-1
                     call MatSetValue(this%A,this%row,ind,val,INSERT_VALUES,ierr)
                     val=-1.*insche%ddamp(i,1)*insche%dloop(3,irun)
                     ind=j+inmesh%blocks_per_param-1
                     call MatSetValue(this%A,this%row,ind,val,INSERT_VALUES,ierr)
                     call VecSetValue(this%b,this%row,rhs,INSERT_VALUES,ierr)
                     call VecSetValue(this%mask_dmp,this%row,rhs,INSERT_VALUES,ierr)
                     call VecSetValue(this%mask_dat,this%row,rhs,INSERT_VALUES,ierr)
                     this%row=this%row+1
                     
                     ! damp the difference between parameters 4 and 3 
                     if ( inopts%npars == 4 ) then
                        val=insche%ddamp(i,2)*insche%dloop(3,irun)
                        ind=j+2*inmesh%blocks_per_param-1 ! zero based in petsc
                        call MatSetValue(this%A,this%row,ind,val,INSERT_VALUES,ierr)
                        val=-1.*insche%ddamp(i,2)*insche%dloop(3,irun)
                        ind=j+3*inmesh%blocks_per_param-1
                        call MatSetValue(this%A,this%row,ind,val,INSERT_VALUES,ierr)
                        call VecSetValue(this%b,this%row,rhs,INSERT_VALUES,ierr)
                        call VecSetValue(this%mask_dmp,this%row,rhs,INSERT_VALUES,ierr)
                        call VecSetValue(this%mask_dat,this%row,rhs,INSERT_VALUES,ierr)
                        this%row=this%row+1
                     end if
                     
                  end do
                  k=k+inmesh%blocks_per_layer(i)
                  
               enddo
            end if
          end subroutine apply_ddamp
    !========================================================



    !========================================================
    !
    !   attach terms related to vp-to-vs scaling
    !
          subroutine apply_scale(this,inopts,inmesh,insche,irun) 
            
            implicit none

            class(matr) :: this
            class(opts) :: inopts
            class(mesh) :: inmesh
            class(sche) :: insche

            PetscScalar   val,rhs
            PetscInt      ind,irun
            PetscInt      i,j,k

            k   = 0
            rhs = 0.d0

            if ( verbosity > 1 ) call PetscPrintf(PETSC_COMM_WORLD,&
                 "Starting with SDAMP, at row "//&
                 trim(int2str(this%row))//"\n",ierr)

            if ( this%processor == 0 ) then ! only estimate scale operator with proc 0              
               do i=1,inopts%nlays
                  if ( verbosity > 1 ) print*,"scale, layer: ",i
                  do j=k+1,k+inmesh%blocks_per_layer(i) ! normal indices
                     select case (inopts%npars)
                     case (4) ! 4 parameters: usually vph,vpv,vsh,vsv
                        ! scale vsh-vph
                        val=insche%dloop(4,irun)*insche%scale(i,3)
                        ind=j+2*inmesh%blocks_per_param-1
                        call MatSetValue(this%A,this%row,ind,val,INSERT_VALUES,ierr)
                        val=-1.*insche%scale(i,1)*insche%dloop(4,irun)*insche%scale(i,3)
                        ind=j-1
                        call MatSetValue(this%A,this%row,ind,val,INSERT_VALUES,ierr)
                        call VecSetValue(this%b,this%row,rhs,INSERT_VALUES,ierr)
                        call VecSetValue(this%mask_dmp,this%row,rhs,INSERT_VALUES,ierr)
                        call VecSetValue(this%mask_dat,this%row,rhs,INSERT_VALUES,ierr)
                        this%row=this%row+1
                        ! scale vsv-vpv
                        val=insche%dloop(4,irun)*insche%scale(i,3)
                        ind=j+3*inmesh%blocks_per_param-1
                        call MatSetValue(this%A,this%row,ind,val,INSERT_VALUES,ierr)
                        val=-1.*insche%scale(i,2)*insche%dloop(4,irun)*insche%scale(i,3) 
                        ind=j+inmesh%blocks_per_param-1
                        call MatSetValue(this%A,this%row,ind,val,INSERT_VALUES,ierr)
                        call VecSetValue(this%b,this%row,rhs,INSERT_VALUES,ierr)
                        call VecSetValue(this%mask_dmp,this%row,rhs,INSERT_VALUES,ierr)
                        call VecSetValue(this%mask_dat,this%row,rhs,INSERT_VALUES,ierr)
                        this%row=this%row+1
                     case (2) ! 2 parameters: usually vp,vs
                        ! scale vp-vs
                        val=insche%dloop(4,irun)*insche%scale(i,3)
                        ind=j+inmesh%blocks_per_param-1
                        call MatSetValue(this%A,this%row,ind,val,INSERT_VALUES,ierr)
                        val=-1.*insche%scale(i,1)*insche%dloop(4,irun)*insche%scale(i,3)
                        ind=j-1
                        call MatSetValue(this%A,this%row,ind,val,INSERT_VALUES,ierr)
                        call VecSetValue(this%b,this%row,rhs,INSERT_VALUES,ierr)
                        call VecSetValue(this%mask_dmp,this%row,rhs,INSERT_VALUES,ierr)
                        call VecSetValue(this%mask_dat,this%row,rhs,INSERT_VALUES,ierr)
                        this%row=this%row+1                   
                     end select
                  end do
                  k=k+inmesh%blocks_per_layer(i)                  
               end do
            end if
          end subroutine apply_scale
    !========================================================




    !========================================================
    !
    !   loops over matrix schedule and assembles A matrix
    !
          subroutine read_submatrices(this,insche) !,Istart,Iend)
            
            implicit none

            class(matr) :: this
            class(sche) :: insche

            ! All data is saved in single precision
            ! but for the global system we need doubles
            real*4,      allocatable :: rhs_in(:)
            real*4,      allocatable :: val_in(:)
            PetscInt,    allocatable :: ind_in(:)
            PetscInt,    allocatable :: pnt_in(:)
            
            PetscScalar, allocatable :: val_petsc(:)
            PetscInt,    allocatable :: ind_petsc(:)
            
            PetscScalar  rhs_petsc
            PetscScalar  cut_weight
            PetscScalar  absttime

            PetscInt     i,j,k,l
            PetscInt     row_total,row
            PetscScalar  valu, weight

            allocate ( pnt_in ( datmx ) )
            allocate ( rhs_in ( datmx ) )
            allocate ( ind_in ( nprmx ) )
            allocate ( val_in ( nprmx ) )
            allocate ( ind_petsc ( nprmx ) )
            allocate ( val_petsc ( nprmx ) )

            this%row = 0
            row_total = 0
            pnt_in(0) = 0

            ! loop over list of datasets
            do j=1,insche%mats_total
            
               open(10,file=trim(insche%path_info(j,3)),status='old')
               open(20,file=trim(insche%path_info(j,4)),status='old')

               ! Read rhs and 
               do i=1,int(insche%meta_info(j,4)) ! loop over pointer values
                  read(10,*) pnt_in(i) !
                  read(20,*) rhs_in(i) ! 
               enddo

               close(10)
               close(20)
               
               insche%fromto_info(j,1) = this%row

               ! Second loop over measurements, insert row by row
               do i=1,int(insche%meta_info(j,4))
               

                  ! Check if current proc owns this row
                  if (((this%row).ge.this%row_start).and.&
                      ((this%row).le.this%row_end-1)) then

                     ! open val and ind vectors
                     open(30,file=trim(insche%path_info(j,1)),status='old',&
                        access='direct',form='unformatted',recl=4)
                     open(40,file=trim(insche%path_info(j,2)),status='old',&
                        access='direct',form='unformatted',recl=4)

                     l=0
                     do k=pnt_in(i-1)+1,pnt_in(i)
                        l=l+1
                        read(30,rec=k) val_in(l) ! 
                        read(40,rec=k) ind_in(l) ! 
                     end do

                     ! estimate weights, from lapos code
                     absttime = abs(rhs_in(i))
                     cut_weight = 1.
                     if ( absttime > insche%meta_info(j,2) ) then
                        cut_weight = exp( insche%meta_info(j,2) - absttime )
                     end if
                     if ( absttime < insche%meta_info(j,3) ) then
                        cut_weight = exp(-insche%meta_info(j,3) + absttime )
                     end if
                  
                     ! Set rhs vector
                     weight = insche%meta_info(j,1) * cut_weight
                     rhs_petsc      = dble ( rhs_in(i) )   * weight
                     val_petsc(1:l) = dble ( val_in(1:l) ) * weight
                     ind_petsc(1:l) = ind_in(1:l) - 1 ! PETSc index starts at 0


                     ! Set one row at a time, in case this row is locally owned
                     call MatSetValues ( this%A, 1, this%row, l, ind_petsc(1:l), &
                                         val_petsc(1:l), INSERT_VALUES, ierr)
                                     
                     ! Set values in parallel vector
                     call VecSetValue  ( this%b, this%row, rhs_petsc, INSERT_VALUES, ierr )
                     call VecSetValue  ( this%mask_dmp, this%row, 1.d0, INSERT_VALUES, ierr )
                     call VecSetValue  ( this%mask_dat, this%row, 0.d0, INSERT_VALUES, ierr )
                     
                     if (verbosity > 1) then
                        if (mod(int(this%row),int(real(this%row_end)/100.)).eq.0) then
                           print*,"Proc ",this%processor, " read", &
                                int((real(this%row)/real(this%row_end))*100)," %"
                        end if
                     end if

                     close(30)
                     close(40)

                  end if ! end if check ownership

                  ! Global row index, PETSc uses 0-based indices
                  this%row = this%row + 1 

               end do

               ! Submatrix 
               insche%fromto_info(j,2) = this%row - 1

               call PetscPrintf(PETSC_COMM_WORLD,"MATRIX: Read submatrix "//&
                                 trim(insche%path_info(j,3))//"\n",ierr)
               !row_total = row_total + int(insche%meta_info(j,4))                    
            end do

            ! free temporary arrays
            deallocate ( pnt_in )
            deallocate ( rhs_in )
            deallocate ( ind_in )
            deallocate ( val_in )
            deallocate ( ind_petsc )
            deallocate ( val_petsc )
            
          end subroutine read_submatrices
    !========================================================


    !========================================================
    !
    !   procedure to assemble a matrix
    !
          subroutine assemble_matrix(this)
            
            implicit none
            class(matr) :: this

            call MatAssemblyBegin ( this%A, MAT_FINAL_ASSEMBLY, ierr)
            call MatAssemblyEnd   ( this%A, MAT_FINAL_ASSEMBLY, ierr)

            call VecAssemblyBegin ( this%b, ierr)
            call VecAssemblyEnd   ( this%b, ierr)

            call VecAssemblyBegin ( this%b_synth, ierr)
            call VecAssemblyEnd   ( this%b_synth, ierr)

            call VecAssemblyBegin ( this%b_adotx, ierr)
            call VecAssemblyEnd   ( this%b_adotx, ierr)

            call VecAssemblyBegin ( this%b_rough, ierr)
            call VecAssemblyEnd   ( this%b_rough, ierr)

            call VecAssemblyBegin ( this%x_synth, ierr)
            call VecAssemblyEnd   ( this%x_synth, ierr)

            call VecAssemblyBegin ( this%mask_dmp, ierr)
            call VecAssemblyEnd   ( this%mask_dmp, ierr)

            call VecAssemblyBegin ( this%mask_dat, ierr)
            call VecAssemblyEnd   ( this%mask_dat, ierr)

            call VecAssemblyBegin ( this%mask_sub, ierr)
            call VecAssemblyEnd   ( this%mask_sub, ierr)

          end subroutine assemble_matrix
    !========================================================


    !========================================================
    !
    !   procedure to assemble all vectors
    !
          subroutine assemble_vectors(this)
            
            implicit none
            class(matr) :: this

            call VecAssemblyBegin ( this%b, ierr)
            call VecAssemblyEnd   ( this%b, ierr)

            call VecAssemblyBegin ( this%b_synth, ierr)
            call VecAssemblyEnd   ( this%b_synth, ierr)

            call VecAssemblyBegin ( this%b_adotx, ierr)
            call VecAssemblyEnd   ( this%b_adotx, ierr)

            call VecAssemblyBegin ( this%b_dummy, ierr)
            call VecAssemblyEnd   ( this%b_dummy, ierr)

            call VecAssemblyBegin ( this%b_rough, ierr)
            call VecAssemblyEnd   ( this%b_rough, ierr)

            call VecAssemblyBegin ( this%x_synth, ierr)
            call VecAssemblyEnd   ( this%x_synth, ierr)

            call VecAssemblyBegin ( this%mask_dmp, ierr)
            call VecAssemblyEnd   ( this%mask_dmp, ierr)

            call VecAssemblyBegin ( this%mask_dat, ierr)
            call VecAssemblyEnd   ( this%mask_dat, ierr)

            call VecAssemblyBegin ( this%mask_sub, ierr)
            call VecAssemblyEnd   ( this%mask_sub, ierr)

          end subroutine assemble_vectors
    !========================================================


    !========================================================
          subroutine read_synth_model(this,inopts,inmesh,insche)

            implicit none

            class(matr) :: this
            class(opts) :: inopts
            class(mesh) :: inmesh
            class(sche) :: insche

            type(mesh) check_mesh ! setup new mesh for checkerboard
            type(opts) check_opts ! setup new opts for checkerboard

            PetscInt,          parameter :: fh=20 ! file handler
            PetscInt                        dummy
            PetscChar(256)                  mode
            PetscInt                        ipar,ios
            PetscInt                        npoints,ipoint
            PetscScalar                     lat,lon,dep,xval
            PetscInt                        i,j,k,h,l,u
            PetscInt                        row_petsc
            PetscBool                       foundlay
            PetscInt                        found

            PetscScalar                     eqinc_check
            PetscScalar,     allocatable :: check_vals(:)

            PetscScalar,     allocatable :: switch_facts(:,:)
            PetscScalar,     allocatable :: valtot(:,:)
            PetscScalar,     allocatable :: numper(:,:)
            
            call PetscPrintf(PETSC_COMM_WORLD,"    reading synthetic input model!\n",ierr)

            if ( this%processor == 0 ) then

               open(fh,file=trim(inopts%synth))               
               read(unit=fh,fmt=*,iostat=ios) mode
               call PetscPrintf(PETSC_COMM_WORLD,"    File: "//trim(inopts%synth)// "\n",ierr)

               select case(trim(mode))

               case ('CHECKERBOARD')
                  
                  !
                  ! Checkerboard model files have the following format
                  !
                  ! CHECKERBOARD # mode in first line
                  ! 20.0         # checkerboard eqincr
                  ! 2.0    -2.0  # Amplitude factors for the n parameters
                  ! 2.0    -2.0  # n = number of columns
                  ! ...
                  !

                  allocate(switch_facts(inopts%nlays,inopts%npars))
                  allocate(valtot(n0max,inopts%nlays))

                  read(unit=fh,fmt=*,iostat=ios) eqinc_check
                  do i=1,inopts%nlays
                     read(fh, fmt=*, iostat=ios) switch_facts(i,1:inopts%npars)
                  end do

                  ! Some more parameters needed to set up checkerboard
                  check_opts%adapt=.false.
                  check_opts%nlays=inopts%nlays
                  check_opts%npars=1 ! not strictly needed
                  check_opts%eqinc=eqinc_check
                  check_opts%eqinc_ref=eqinc_check ! this is important

                  ! Setup new checkerboard mesh
                  call PetscPrintf(PETSC_COMM_WORLD,"...\n",ierr)
                  call check_mesh%setup_mesh(check_opts)
                  call PetscPrintf(PETSC_COMM_WORLD,"...\n",ierr)

                  ! Initialize checkerboard values
                  allocate(check_vals(check_mesh%blocks_per_layer(1)))
                  check_vals(1)=1.d0
                  do i=2,check_mesh%blocks_per_layer(1) ! just do it for layer 1
                     check_vals(i) = check_vals(i-1)*(-1.d0)
                     if(check_mesh%xlamin(i,1).lt.check_mesh%xlamin(i-1,1)) then                        
                        check_vals(i) = check_vals(i)*(-1.d0)
                     end if
                  end do
                    
                  ! Loop over parameter
                  do u=1,inopts%npars
                     do i=1,inopts%nlays
                        do j=1,inmesh%blocks_per_layer(i)
                           found=0
                           do h=1,check_mesh%blocks_per_layer(1)
                              if ((inmesh%locent(j,i).ge.check_mesh%xlomin(h,1)).and.&
                                   (inmesh%locent(j,i).lt.check_mesh%xlomax(h,1)).and.&
                                   (inmesh%lacent(j,i).ge.check_mesh%xlamin(h,1)).and.&
                                   (inmesh%lacent(j,i).lt.check_mesh%xlamax(h,1))) then
                                 valtot(j,i)=check_vals(h)*switch_facts(i,u)
                                 found=found+1
                              end if
                           end do

                           ! Some sanity checks
                           if (found.lt.1.) then
                              call PetscPrintf(PETSC_COMM_WORLD,"Error in checkerboard creation",ierr)
                              stop
                           else if (found.gt.1) then
                              call PetscPrintf(PETSC_COMM_WORLD,"Error in checkerboard creation",ierr)
                              stop
                           end if
                        end do
                     end do ! end fo loop over inopts%nlay

                     ! Set values of x_synth
                     row_petsc=(u-1)*inmesh%blocks_per_param
                     do j=1,inmesh%nlays                        
                        do k=1,inmesh%blocks_per_layer(j)
                           ! print*,valtot(k,j)
                           call VecSetValue(this%x_synth,row_petsc,&
                                valtot(k,j)/100.d0,INSERT_VALUES,ierr)
                           row_petsc=row_petsc+1
                        enddo
                     enddo

                  end do ! end of loop over npars
                                                                          
                  deallocate(switch_facts)                 
                  deallocate(valtot)

               case ('POINTCLOUD')

                  !
                  ! Pointcloud model files have the following format
                  !
                  ! POINTCLOUD # mode in first line
                  ! 1000000    # number of points per parameter
                  !  1.275     # npar*npoints model coefficinets in %
                  !  1.556
                  !  0.988
                  ! -0.899
                  ! ...
                  !

                  allocate(numper(n0max,inopts%nlays))
                  allocate(valtot(n0max,inopts%nlays))
                  read(unit=fh,fmt=*,iostat=ios) npoints

                  ! Loop over parameter
                  do u=1,inopts%npars

                     ! Initialize arrays for new parameter
                     numper=0.d0
                     valtot=0.d0
                     l=0
                     
                     ! Loop over npoints per parameter
                     do ipoint=1,npoints

                        ! Read point
                        read(unit=fh,fmt=*,iostat=ios) lat,lon,dep,xval

                        ! Shift grid
                        lon=lon+360.d0 
                        if (lon.gt.(359.60)) lon=lon-360.d0 

                        ! In which layer are we
                        foundlay=.false.

                        do h=1,inmesh%nlays
                           if((dep.gt.insche%layer(h)).and.&
                                (dep.le.insche%layer(h+1)))then
                              if(.not.foundlay)then
                                 foundlay=.true.
                                 l=h
                              else
                                 call PetscPrintf(PETSC_COMM_WORLD,"Error in read_synth_model",ierr)
                                 call PetscPrintf(PETSC_COMM_WORLD,"Point lies b/w two layers",ierr)
                                 stop 
                              endif
                           endif
                        enddo
                        if (l==0) stop "Something went wrong"
                                               
                        ! In which pixel are we
                        do i=1,inmesh%blocks_per_layer(l)
                           if ((lon.ge.inmesh%xlomin(i,l)).and.&
                                (lon.le.inmesh%xlomax(i,l)).and.&
                                (lat.ge.inmesh%xlamin(i,l)).and.&
                                (lat.le.inmesh%xlamax(i,l))) then                           

                              numper(i,l)=numper(i,l)+1.d0
                              valtot(i,l)=valtot(i,l)+xval

                           end if
                        end do
                     end do ! end of loop over npoints

                     ! Set values of x_synth
                     row_petsc=(u-1)*inmesh%blocks_per_param
                     do j=1,inmesh%nlays
                        do k=1,inmesh%blocks_per_layer(j)
                           valtot(k,j)=valtot(k,j)/numper(k,j)
                           ! print*,valtot(k,j)
                           call VecSetValue(this%x_synth,row_petsc,&
                                valtot(k,j)/100.d0,INSERT_VALUES,ierr)
                           row_petsc=row_petsc+1
                        enddo
                     enddo

                  end do ! end of loop over npar

                  deallocate(numper)
                  deallocate(valtot)
                  
               case ('VOXEL')

                  !
                  ! Voxel model files have the following format
                  !
                  ! VOXEL     # mode in first line
                  ! 1  1.275  # n voxel id (int) and coeff in %
                  ! 2  1.556  # n = nvoxels * nparams
                  ! 3  0.988 
                  ! 4 -0.899
                  ! ...
                  !

                  do ipar=1,inmesh%blocks_all_param
                     read(unit=fh,fmt=*,iostat=ios) dummy,xval            
                     call VecSetValue(this%x_synth,ipar-1,&
                          xval/100.d0,INSERT_VALUES,ierr)
                  end do                   
               end select

               close(fh)  
            end if
            
          end subroutine read_synth_model
    !========================================================


    !========================================================
          subroutine store_rhs(this)
            implicit none
            class(matr) :: this
            call PetscPrintf(PETSC_COMM_WORLD,"    storing rhs vector!\n",ierr)
            call VecCopy(this%b,this%b_store,ierr)
          end subroutine store_rhs
    !========================================================


    !========================================================
          subroutine restore_rhs(this)
            implicit none
            class(matr) :: this
            call PetscPrintf(PETSC_COMM_WORLD,"    re-storing rhs vector!\n",ierr)
            call VecCopy(this%b_store,this%b,ierr)
          end subroutine restore_rhs
    !========================================================


    !========================================================
          subroutine compute_synthetics(this)
            implicit none
            class(matr) :: this
            call PetscPrintf(PETSC_COMM_WORLD,"    computing synthetics!\n",ierr)
            call MatMult(this%A,this%x_synth,this%b_synth,ierr)           
            call VecPointwiseMult(this%b,this%b_synth,this%mask_dmp,ierr)            
          end subroutine compute_synthetics
    !========================================================


    !========================================================
          subroutine compute_adotx(this,xvec,bvec,mask)

            implicit none
            class(matr)        :: this

            Vec, intent(inout) :: xvec
            Vec, intent(inout) :: bvec

            Vec, intent(in), optional :: mask

            call PetscPrintf(PETSC_COMM_WORLD,"    computing A dot x!\n",ierr)            
            call MatMult(this%A,xvec,bvec,ierr)

            if (present(mask)) then
               call VecPointwiseMult(bvec,bvec,mask,ierr)
            end if
            
          end subroutine compute_adotx
    !========================================================


    !========================================================
          subroutine destroy_matrix(this)
            implicit none
            class(matr) :: this
            call VecDestroy(this%x,ierr)
            call VecDestroy(this%b,ierr)
            call MatDestroy(this%A,ierr)           
          end subroutine destroy_matrix
    !========================================================



      end module module_matrix

!============================================================




!============================================================
! 
!  module containing everything related to the solver
!
      module module_solver

          use global_param     
          use module_options
          use module_matrix

          implicit none

          type,public :: solv
             KSP            ksp ! kralov subspace solver context  
             PC             pc  ! preconditioner context
           contains
             procedure :: solve_system
             procedure :: destroy_solver
          end type solv

        contains


    !========================================================
    !
    !   solve the system
    !
          subroutine solve_system(this,inmatr,inopts)
            
            implicit none
            class(solv) :: this
            class(matr) :: inmatr
            class(opts) :: inopts
 
            PetscInt   mxit
            PetscReal  atol
            PetscReal  dtol
            PetscReal  rtol

            ! GMRES ON THE NORMAL EQUATIONS
            if ( inopts%type == 'normal' .or. &
                 inopts%type == 'gmres' ) then

               ! Setup Normal equations
               call MatCreateNormal ( inmatr%A, inmatr%ATA, ierr )
               call VecDuplicate ( inmatr%x, inmatr%ATb, ierr )
               call MatMultTranspose ( inmatr%A, inmatr%b, inmatr%ATb, ierr )

               ! Setup gmres solver context
               call KSPCreate ( PETSC_COMM_WORLD, this%ksp, ierr )
               call KSPSetOperators ( this%ksp, inmatr%ATA, inmatr%ATA, &
                                      SAME_NONZERO_PATTERN, ierr )

               ! Default KSP Type is 
               call KSPSetType ( this%ksp, KSPGMRES, ierr )
               call MatSetUp ( inmatr%ATA, ierr )

               ! Default Preconditioner
               call KSPGetPC ( this%ksp, this%pc, ierr)
               call PCSetType  (this%pc, PCNONE, ierr)

               ! Allow to set KSP options at runtime
               call KSPSetFromOptions ( this%ksp, ierr)               

!               call PetscPrintf(PETSC_COMM_WORLD,"SOLVER: Setup a KSP context for "//&
!                    //"a\n",ierr)

               ! Solve system
               call KSPSolve ( this%ksp, inmatr%ATb, inmatr%x, ierr)

               
               

            ! LSQR ON THE RECTANGULAR SYSTEM
            else if ( inopts%type == 'rect' .or. &
                      inopts%type == 'lsqr' ) then

               ! Setup solver context
               call KSPCreate ( PETSC_COMM_WORLD, this%ksp, ierr )

               !  Set operators. Here the matrix that defines the linear system
               !  also serves as the preconditioning matrix.
               call KSPSetOperators ( this%ksp, inmatr%A, inmatr%A, &
                                      DIFFERENT_NONZERO_PATTERN, ierr )

               ! Set KSP Type
               call KSPSetType ( this%ksp, KSPLSQR, ierr )

               ! Set preconditioner
               call KSPGetPC   ( this%ksp, this%pc ,ierr )
               call PCSetType  ( this%pc, PCNONE, ierr )

               ! Allow to set KSP options at runtime
               call KSPSetFromOptions ( this%ksp, ierr )

               ! Solve system
               call KSPSolve (this%ksp, inmatr%b, inmatr%x, ierr)

            end if
            
          end subroutine solve_system
    !========================================================



    !========================================================
          subroutine destroy_solver(this)
            implicit none
            class(solv) :: this
            call KSPDestroy(this%ksp,ierr)
          end subroutine destroy_solver
    !========================================================



      end module module_solver

!============================================================




!============================================================
!  
!   module associated with background models
!
        module module_bgmod

          use global_param     

          implicit none

          type,public :: back
             PetscScalar      rho
             PetscScalar      vph
             PetscScalar      vpv
             PetscScalar      vsh
             PetscScalar      vsv
             PetscScalar      qka
             PetscScalar      qmu
           contains
             procedure, pass :: get_model_coeff
             procedure, pass :: get_model_coeff_aver
             procedure, nopass :: get_idom
             procedure, nopass :: prem_ani
          end type back

        contains

    !========================================================
    !
    !  function returning average model coeff for radius interval
    !
          doubleprecision function get_model_coeff_aver(this,dbeg,dend,incr,param,bgmod)

            class(back),      intent(in) :: this
            PetscScalar,      intent(in) :: dbeg
            PetscScalar,      intent(in) :: dend
            PetscScalar,      intent(in) :: incr
            character(len=*), intent(in) :: bgmod
            character(len=*), intent(in) :: param

            PetscScalar radius
            PetscScalar coeff
            PetscScalar count
            PetscScalar rbeg
            PetscScalar rend
            PetscInt idom

            rbeg = r_earth - dbeg 
            rend = r_earth - dend
            
            coeff = 0.d0
            count = 0.d0

            radius = rbeg
            do while (radius > rend)
               radius = radius - incr
               count  = count + 1.d0
               idom   = this%get_idom(radius,bgmod) 
               select case(trim(bgmod))
               case('prem_iso')             
                  coeff = coeff + prem_ani(radius,param,idom)
               case('prem_ani')
                  coeff = coeff + prem_ani(radius,param,idom)
               case default
                  write(6,*) 'POSTPROC: Unknown background model: ', trim(bgmod)
                  stop
               end select  
            end do

            get_model_coeff_aver = coeff / count

          end function get_model_coeff_aver
    !========================================================


    !========================================================
    !
    !  function returning model coefficients at specified depths
    !
          doubleprecision function get_model_coeff(this,depth,param,bgmod)

            class(back),      intent(in) :: this
            PetscScalar,      intent(in) :: depth
            character(len=*), intent(in) :: bgmod
            character(len=*), intent(in) :: param
            PetscScalar radius
            PetscInt idom

            radius = r_earth - depth
            idom = this%get_idom(radius,bgmod)            
            select case(trim(bgmod))
            case('prem_iso')             
               get_model_coeff = prem_ani(radius,param,idom) !@TODO prem_iso              
            case('prem_ani')
               get_model_coeff = prem_ani(radius,param,idom)
            case default
               write(6,*) 'POSTPROC: Unknown background model: ', trim(bgmod)
               stop
             end select

          end function get_model_coeff
    !========================================================


    !========================================================
          doubleprecision function prem_ani(radius,param,idom)           

            PetscScalar,  intent(in) :: radius
            character(len=*), intent(in) :: param
            PetscInt,     intent(in) :: idom

            PetscScalar vpv_prem, vph_prem
            PetscScalar vsh_prem, vsv_prem 
            PetscScalar ro_prem, eta_aniso
            PetscScalar qmu, qkappa
            PetscScalar x_prem
            
            ! Initialize variables
            x_prem = radius / r_earth ! Radius (normalized to x(surface)=1)
            eta_aniso = 1.d0
            ro_prem   = 0.d0
            vpv_prem  = 0.d0
            vph_prem  = 0.d0
            vsv_prem  = 0.d0
            vsh_prem  = 0.d0
            qmu       = 0.d0
            qkappa    = 0.d0
           
            select case (idom)
            case (1)
               ro_prem  = 2.6
               vpv_prem = 5.8
               vsv_prem = 3.2
               vph_prem = vpv_prem
               vsh_prem = vsv_prem
               qmu = 600.0
               qkappa = 57827.0
            case (2) ! lower crustal layer
               ro_prem  = 2.9
               vpv_prem = 6.8
               vsv_prem = 3.9
               vph_prem = vpv_prem
               vsh_prem = vsv_prem
               qmu = 600.0
               qkappa = 57827.0
            case (3) ! upper mantle
               ro_prem   =  2.6910 + 0.6924 * x_prem
               vpv_prem  =  0.8317 + 7.2180 * x_prem
               vph_prem  =  3.5908 + 4.6172 * x_prem
               vsv_prem  =  5.8582 - 1.4678 * x_prem
               vsh_prem  = -1.0839 + 5.7176 * x_prem
               eta_aniso =  3.3687 - 2.4778 * x_prem
               qmu = 600.0
               qkappa = 57827.0
            case (4)
               ro_prem  =  7.1089 -  3.8045 * x_prem
               vpv_prem = 20.3926 - 12.2569 * x_prem
               vsv_prem =  8.9496 -  4.4597 * x_prem
               vph_prem = vpv_prem
               vsh_prem = vsv_prem
               qmu = 143.0
               qkappa = 57827.0
            case (5)
               ro_prem  = 11.2494 -  8.0298 * x_prem
               vpv_prem = 39.7027 - 32.6166 * x_prem
               vsv_prem = 22.3512 - 18.5856 * x_prem
               vph_prem = vpv_prem
               vsh_prem = vsv_prem
               qmu = 143.0
               qkappa = 57827.0
            case (6)
               ro_prem  =  5.3197 - 1.4836 * x_prem
               vpv_prem = 19.0957 - 9.8672 * x_prem
               vsv_prem =  9.9839 - 4.9324 * x_prem
               vph_prem = vpv_prem
               vsh_prem = vsv_prem
               qmu = 143.0
               qkappa = 57827.0
            case (7) ! lower mantle
               ro_prem  =  7.9565 - 6.4761 * x_prem + &
                    5.5283 * x_prem**2 - 3.0807 * x_prem**3
               vpv_prem = 29.2766 -23.6027 * x_prem + &
                    5.5242 * x_prem**2 - 2.5514 * x_prem**3
               vsv_prem = 22.3459 -17.2473 * x_prem - &
                    2.0834 * x_prem**2 + 0.9783 * x_prem**3
               vph_prem = vpv_prem
               vsh_prem = vsv_prem
               qmu = 312.0
               qkappa = 57827.0
            case (8)
               ro_prem  =  7.9565 -  6.4761 * x_prem + &
                    5.5283 * x_prem**2 -  3.0807 * x_prem**3
               vpv_prem = 24.9520 - 40.4673 * x_prem + &
                    51.4832 * x_prem**2 - 26.6419 * x_prem**3
               vsv_prem = 11.1671 - 13.7818 * x_prem + &
                    17.4575 * x_prem**2 -  9.2777 * x_prem**3
               vph_prem = vpv_prem
               vsh_prem = vsv_prem
               qmu = 312.0
               qkappa = 57827.0
            case (9)
               ro_prem  =  7.9565 - 6.4761 * x_prem + &
                    5.5283 * x_prem**2 - 3.0807 * x_prem**3
               vpv_prem = 15.3891 - 5.3181 * x_prem + &
                    5.5242 * x_prem**2 - 2.5514 * x_prem**3
               vsv_prem =  6.9254 + 1.4672 * x_prem - &
                    2.0834 * x_prem**2 + 0.9783 * x_prem**3
               vph_prem = vpv_prem
               vsh_prem = vsv_prem
               qmu = 312.0
               qkappa = 57827.0
            case (10)
               ro_prem  = 12.5815 - 1.2638 * x_prem - &
                    3.6426 * x_prem**2 -  5.5281 * x_prem**3
               vpv_prem = 11.0487 - 4.0362 * x_prem + &
                    4.8023 * x_prem**2 - 13.5732 * x_prem**3
               vsv_prem =  0.0
               vph_prem = vpv_prem
               vsh_prem = vsv_prem
               qmu = 0.0
               qkappa = 57827.0
            case (11)
               ro_prem  = 13.0885 - 8.8381 * x_prem**2
               vpv_prem = 11.2622 - 6.3640 * x_prem**2
               vsv_prem =  3.6678 - 4.4475 * x_prem**2
               vph_prem = vpv_prem
               vsh_prem = vsv_prem
               qmu = 84.6
               qkappa = 1327.7
            end select
            
            select case (trim(param))
            case ('rho')
               prem_ani = ro_prem * 1000.
            case ('vpv')
               prem_ani = vpv_prem * 1000.
            case ('vsv')
               prem_ani = vsv_prem * 1000.
            case ('vph')
               prem_ani = vph_prem * 1000.
            case ('vsh')
               prem_ani = vsh_prem * 1000.
            case ('eta')
               prem_ani = eta_aniso
            case ('qmu')
               prem_ani = qmu
            case ('qka')
               prem_ani = qkappa
            case default 
               call PetscPrintf(PETSC_COMM_WORLD,"POSTPROC: Error in prem_ani",ierr)
               stop
            end select

          end function prem_ani
    !========================================================


    !========================================================
          integer function get_idom(r0,bgmod)           
            
            character(len=*), intent(in) :: bgmod
            PetscScalar, intent(in)      :: r0
            PetscScalar discont(30)
            PetscInt idom
            PetscInt ndisc

            select case (trim(bgmod))
            case ('prem_ani')
               ndisc = 11
               discont(1) = 6371.    ! UPPER CRUST
               discont(2) = 6356.    ! LOWER CRUST
               discont(3) = 6346.6   ! MOHO --> 220
               discont(4) = 6151.    ! TRANSITION ZONE: 220 --> 400
               discont(5) = 5971.    ! 400 --> 600
               discont(6) = 5771.    ! 600 --> 670 
               discont(7) = 5701.    ! 670 --> 770
               discont(8) = 5600.    ! LOWER MANTLE: 770 --> TOP D"
               discont(9) = 3630.    ! D" LAYER 
               discont(10) = 3480.   ! FLUID OUTER CORE: CMB --> ICB
               discont(11) = 1221.5  ! SOLID INNER CORE: ICB --> CENTER 
               discont(12) = 0.      ! CENTER OF THE EARTH
            case ('prem_iso')
               call PetscPrintf(PETSC_COMM_WORLD,"POSTPROC: Error in get_idom",ierr)
               stop
            case default
               call PetscPrintf(PETSC_COMM_WORLD,"POSTPROC: Error in prem_ani",ierr)
               stop
            end select

            get_idom = 0
            do idom=1,ndisc
               if ((r0.le.discont(idom)).and.&
                   (r0.gt.discont(idom+1))) then
                  get_idom = idom
                  exit
               end if
            end do
            
          end function get_idom

        end module module_bgmod
!============================================================        




!============================================================
!
!  postprocessing: reparameterize and output to disk
!  
      module module_postproc

          use global_param     
          use module_options
          use module_solver
          use module_matrix
          use module_bgmod

          implicit none

          type,public :: post
             Vec              xout

             Vec              rough
             Vec              sol_postproc

             VecScatter       ctx
             Mat              D0
             PetscViewer      viewer
             IS               isrow
             IS               iscol

             ! for reparameterization
             PetscScalar, allocatable :: sol_raw(:)
             PetscScalar, allocatable :: sol_dlnvsh(:)
             PetscScalar, allocatable :: sol_dlnvsv(:)
             PetscScalar, allocatable :: sol_dlnvph(:)
             PetscScalar, allocatable :: sol_dlnvpv(:)
             PetscScalar, allocatable :: sol_dlnvp(:)
             PetscScalar, allocatable :: sol_dlnvs(:)
             PetscScalar, allocatable :: sol_dlnvc(:)
             PetscScalar, allocatable :: sol_vpv(:)
             PetscScalar, allocatable :: sol_vph(:)
             PetscScalar, allocatable :: sol_vsv(:)
             PetscScalar, allocatable :: sol_vsh(:)
             PetscScalar, allocatable :: sol_vp(:)
             PetscScalar, allocatable :: sol_vs(:)
             PetscScalar, allocatable :: sol_xi(:)
             PetscScalar, allocatable :: sol_ph(:)

             ! roughness and model norm
             PetscScalar, allocatable :: rough_abs(:)
             PetscScalar, allocatable :: rough_nrm(:)
             PetscScalar, allocatable :: norm_rms(:)
             PetscScalar, allocatable :: norm_abs(:)
             PetscScalar, allocatable :: varr_cumm(:)
             PetscScalar, allocatable :: varr(:,:)
             PetscInt,    allocatable :: its(:)

             ! for xdmf output
             PetscInt,    allocatable :: connectivity(:,:)
             PetscScalar, allocatable :: vertices(:,:)

             ! booleans for output
             PetscBool do_vsh,do_vsv
             PetscBool do_vph,do_vpv
             PetscBool do_xi,do_ph
             PetscBool do_vs,do_vp
             PetscBool do_vc

             PetscBool initialized


           contains

             procedure, pass :: initialize_postproc
             procedure, pass :: reparam_solution
             procedure, pass :: destroy_postproc
             procedure, pass :: export_solution
             procedure, pass :: save_iterations
             procedure, pass :: compute_norm
             procedure, pass :: compute_rough
             procedure, pass :: compute_varr
             procedure, nopass :: sph2cart
             procedure, nopass :: dump_model_ascii
             procedure, pass :: dump_model_xdmf
             procedure, pass :: dump_run_info
             procedure, pass :: dump_varr

          end type post

        contains



    !========================================================
    !
    !   collect solution vector from all processors
    !
          subroutine initialize_postproc(this,inopts,inmatr,inmesh,insche,sol_vec)

            implicit none
            class(post) :: this
            class(opts) :: inopts
            class(matr) :: inmatr
            class(mesh) :: inmesh
            class(sche) :: insche

            Vec, intent(in), optional :: sol_vec
            PetscInt i

            ! In case a solution vector is passed to init, replace inmatr%x with
            ! this solution vector! Relevant for exporting solutions
            if (present(sol_vec)) then
               call inmatr%assemble_vectors()
               call VecCopy(sol_vec,inmatr%x,ierr)
            end if

            ! collect solution vector from all processors and stores in one sequential vector
            call VecScatterCreateToZero(inmatr%x,this%ctx,this%xout,ierr)
            call VecScatterBegin(this%ctx,inmatr%x,this%xout,INSERT_VALUES,SCATTER_FORWARD,ierr)
            call VecScatterEnd(this%ctx,inmatr%x,this%xout,INSERT_VALUES,SCATTER_FORWARD,ierr)

            call PetscPrintf(PETSC_COMM_WORLD,"    initialize postprocessing!\n",ierr)            
            if (.not. this%initialized ) then

               ! allocate potential output vectors (not all of them may be used)
               if ( inmatr%processor == 0 ) then

                  call PetscPrintf(PETSC_COMM_WORLD,"    allocate buffers (only done once)!\n",ierr)
                  ! Allocate output arrays
                  allocate ( this%sol_raw(inmesh%blocks_all_param) )
                  allocate ( this%sol_dlnvs(inmesh%blocks_per_param) )
                  allocate ( this%sol_dlnvc(inmesh%blocks_per_param) )
                  allocate ( this%sol_dlnvp(inmesh%blocks_per_param) )
                  allocate ( this%sol_dlnvsh(inmesh%blocks_per_param) )
                  allocate ( this%sol_dlnvsv(inmesh%blocks_per_param) )
                  allocate ( this%sol_dlnvph(inmesh%blocks_per_param) )
                  allocate ( this%sol_dlnvpv(inmesh%blocks_per_param) )
                  allocate ( this%sol_vph(inmesh%blocks_per_param) )
                  allocate ( this%sol_vpv(inmesh%blocks_per_param) )
                  allocate ( this%sol_vsh(inmesh%blocks_per_param) )
                  allocate ( this%sol_vsv(inmesh%blocks_per_param) )
                  allocate ( this%sol_vp(inmesh%blocks_per_param) )
                  allocate ( this%sol_vs(inmesh%blocks_per_param) )
                  allocate ( this%sol_xi(inmesh%blocks_per_param) )
                  allocate ( this%sol_ph(inmesh%blocks_per_param) )

                  if ( inopts%format == 'xdmf' ) then
                     allocate ( this%connectivity(8,inmesh%blocks_per_param) ) 
                     allocate ( this%vertices(3,8*inmesh%blocks_per_param) )
                  end if

               end if

               ! Buffers for model roughness, norm and iterations
               ! need to be allocated on all processors
               allocate ( this%its(insche%loop_total) )
               allocate ( this%norm_abs(insche%loop_total) )
               allocate ( this%norm_rms(insche%loop_total) )
               allocate ( this%rough_abs(insche%loop_total) )
               allocate ( this%rough_nrm(insche%loop_total) )
               allocate ( this%varr_cumm(insche%loop_total) )
               allocate ( this%varr(insche%mats_total,insche%loop_total) )

               ! Postprocessing instance is now initialized
               this%initialized = .true.
                              
            end if

            ! This second initialization part is 
            ! repeated after each inversion run
            if ( inmatr%processor == 0 ) then
               ! @TODO, I didn't find a better way to collect my
               ! result from the PETSc vector object, this can be
               ! done more elegantly, for sure
               do i=1,inmesh%blocks_all_param                 
                  call VecGetValues(this%xout,1,i-1,this%sol_raw(i),ierr)
               end do
                  
               ! Initialize output arrays
               this%sol_dlnvs = 0.d0
               this%sol_dlnvc = 0.d0
               this%sol_dlnvp = 0.d0
               this%sol_dlnvsh = 0.d0
               this%sol_dlnvph = 0.d0
               this%sol_dlnvsv = 0.d0
               this%sol_dlnvpv = 0.d0
               this%sol_vp = 0.d0
               this%sol_vs = 0.d0
               this%sol_vph = 0.d0
               this%sol_vpv = 0.d0
               this%sol_vsh = 0.d0
               this%sol_vsv = 0.d0
               this%sol_xi = 0.d0
               this%sol_ph = 0.d0

               ! Initialize logicals
               this%do_vsh = .false.
               this%do_vsv = .false.
               this%do_vph = .false.
               this%do_vpv = .false.
               this%do_vs  = .false.
               this%do_vp  = .false.
               this%do_vc  = .false.
            end if

            call PetscPrintf(PETSC_COMM_WORLD,"    done with initialization!\n",ierr)
           
          end subroutine initialize_postproc
    !========================================================



    !========================================================          
    !
    !   destroys all temporary buffers needed to export model
    !
          subroutine destroy_postproc(this)

            implicit none
            class(post) :: this

               ! Allocate output arrays
               if ( allocated ( this%sol_raw ) )    deallocate ( this%sol_raw )
               if ( allocated ( this%sol_dlnvs ) )  deallocate ( this%sol_dlnvs )
               if ( allocated ( this%sol_dlnvp ) )  deallocate ( this%sol_dlnvp )
               if ( allocated ( this%sol_dlnvsh ) ) deallocate ( this%sol_dlnvsh )
               if ( allocated ( this%sol_dlnvsv ) ) deallocate ( this%sol_dlnvsv )
               if ( allocated ( this%sol_dlnvph ) ) deallocate ( this%sol_dlnvph )
               if ( allocated ( this%sol_dlnvpv ) ) deallocate ( this%sol_dlnvpv )
               if ( allocated ( this%sol_vph ) )    deallocate ( this%sol_vph )
               if ( allocated ( this%sol_vpv ) )    deallocate ( this%sol_vpv )
               if ( allocated ( this%sol_vsh ) )    deallocate ( this%sol_vsh )
               if ( allocated ( this%sol_vsv ) )    deallocate ( this%sol_vsv )
               if ( allocated ( this%sol_vp ) )     deallocate ( this%sol_vp )
               if ( allocated ( this%sol_vs ) )     deallocate ( this%sol_vs )
               if ( allocated ( this%sol_xi ) )     deallocate ( this%sol_xi )
               if ( allocated ( this%sol_ph ) )     deallocate ( this%sol_ph )            
               if ( allocated ( this%connectivity ) ) deallocate ( this%connectivity )
               if ( allocated ( this%vertices ) ) deallocate ( this%vertices )
               if ( allocated ( this%its ) )     deallocate ( this%its )
               if ( allocated ( this%norm_abs ) )     deallocate ( this%norm_abs )
               if ( allocated ( this%norm_rms ) )     deallocate ( this%norm_rms )            
               if ( allocated ( this%rough_nrm ) )     deallocate ( this%rough_nrm )            
               if ( allocated ( this%rough_abs ) )     deallocate ( this%rough_abs )
               if ( allocated ( this%varr_cumm ) )     deallocate ( this%varr_cumm )

               call PetscViewerDestroy(this%viewer,ierr)
               call VecScatterDestroy(this%ctx,ierr)
               call VecDestroy(this%xout,ierr)
               call MatDestroy(this%D0,ierr)           

          end subroutine destroy_postproc
    !========================================================
          


    !========================================================
    !
    !  reparameterize from original parameterization to
    !
          subroutine reparam_solution(this,inopts,inmatr,inmesh,insche)

            implicit none
            
            class(post) :: this          
            class(matr) :: inmatr
            class(opts) :: inopts
            class(sche) :: insche
            class(mesh) :: inmesh

            type(back)  :: my_back ! background model type
            
            PetscScalar dlnvsh,dlnvsv
            PetscScalar dlnvph,dlnvpv
            PetscScalar dlnvp,dlnvs,dlnvc
            PetscScalar vsh,vsv
            PetscScalar vph,vpv
            PetscScalar gamma
            PetscScalar ref_vsh,ref_vsv
            PetscScalar ref_vph,ref_vpv
            PetscScalar ref_vs,ref_vp
            PetscScalar, parameter :: incr = 1.d0 
            PetscInt i,j,k,l,ind
            PetscInt ix_vsv,ix_vsh
            PetscInt ix_vpv,ix_vph

            ix_vsv = 0
            ix_vsh = 0
            ix_vph = 0
            ix_vpv = 0

            call PetscPrintf(PETSC_COMM_WORLD,"    reparameterizing model!\n",ierr)  
            if (inmatr%processor == 0) then              
               do i=1,inopts%npars
                  if ((inopts%parameter(i) == 'vsh').or.&
                      (inopts%parameter(i) == 'VSH')) then
                     if (verbosity > 1) print*,"V_SH solution available"
                     this%do_vsh = .true.
                     ix_vsh = i
                  end if
                  if ((inopts%parameter(i) == 'vsv').or.&
                      (inopts%parameter(i) == 'VSV')) then
                     if (verbosity > 1) print*,"V_SV solution available"
                     this%do_vsv = .true.
                     ix_vsv = i
                  end if
                  if ((inopts%parameter(i) == 'vph').or.&
                      (inopts%parameter(i) == 'VPH')) then
                     if (verbosity > 1) print*,"V_PH solution available"
                     this%do_vph = .true.
                     ix_vph = i
                  end if
                  if ((inopts%parameter(i) == 'vpv').or.&
                      (inopts%parameter(i) == 'VPV')) then
                     if (verbosity > 1) print*,"V_PV solution available"
                     this%do_vpv = .true.
                     ix_vpv = i
                  end if
               end do
               if (this%do_vsv.and.this%do_vsh) then 
                  this%do_vs = .true.
                  this%do_xi = .true.
               end if
               if (this%do_vpv.and.this%do_vph) then 
                  this%do_vp = .true.
                  this%do_ph = .true.
               end if
               if (this%do_vp.and.this%do_vs) then 
                  this%do_vc = .true.
               end if


               k = 0
               do l=2,inopts%nlays+1
                  i = l-1
                  if (verbosity > 1) print*,"Reparameterizing at layer: ",i

                  ! Evaluate arithmetic average of background model within current layer
                  ref_vsh = my_back%get_model_coeff_aver(insche%layer(l-1),insche%layer(l),&
                                                         incr,'vsh',inopts%refmod)
                  ref_vsv = my_back%get_model_coeff_aver(insche%layer(l-1),insche%layer(l),&
                                                         incr,'vsv',inopts%refmod)
                  ref_vph = my_back%get_model_coeff_aver(insche%layer(l-1),insche%layer(l),&
                                                         incr,'vph',inopts%refmod)
                  ref_vpv = my_back%get_model_coeff_aver(insche%layer(l-1),insche%layer(l),&
                                                         incr,'vpv',inopts%refmod)

                  ! Voigt average velocities
                  ref_vp = dsqrt ( ( ref_vpv**2 + 4*ref_vph**2 ) / 5 )
                  ref_vs = dsqrt ( ( 2*ref_vsv**2 + ref_vsh**2 ) / 3 )
                  gamma  = (4.d0/3.d0) * ref_vs**2 / ref_vp**2

                  do j=k+1,k+inmesh%blocks_per_layer(i) 
                     if (this%do_vph) then
                        ind = j + (ix_vph-1)*inmesh%blocks_per_param 
                        this%sol_vph(j) = ref_vph + ref_vph*this%sol_raw(ind)
                        this%sol_dlnvph(j) = this%sol_raw(ind)
                     end if
                     if (this%do_vpv) then
                        ind = j + (ix_vpv-1)*inmesh%blocks_per_param 
                        this%sol_vpv(j) = ref_vpv + ref_vpv*this%sol_raw(ind)
                        this%sol_dlnvpv(j) = this%sol_raw(ind)
                     end if
                     if (this%do_vsh) then
                        ind = j + (ix_vsh-1)*inmesh%blocks_per_param 
                        this%sol_vsh(j) = ref_vsh + ref_vsh*this%sol_raw(ind)
                        this%sol_dlnvsh(j) = this%sol_raw(ind)
                     end if                     
                     if (this%do_vsv) then
                        ind = j + (ix_vsv-1)*inmesh%blocks_per_param 
                        this%sol_vsv(j) = ref_vsv + ref_vsv*this%sol_raw(ind) 
                        this%sol_dlnvsv(j) = this%sol_raw(ind)
                     end if
                     if (this%do_vp) then
                        this%sol_vp(j) = dsqrt ( ( this%sol_vpv(j)**2 + 4*this%sol_vph(j)**2 ) / 5 )
                        this%sol_dlnvp(j) = ( this%sol_vp(j) - ref_vp ) / ref_vp
                     end if
                     if (this%do_vs) then
                        this%sol_vs(j) = dsqrt ( ( 2*this%sol_vsv(j)**2 + this%sol_vsh(j)**2 ) / 3 )
                        this%sol_dlnvs(j) = ( this%sol_vs(j) - ref_vs ) / ref_vs
                     end if
                     if (this%do_vc) then
                        ! This equation is taken from Paulas paper
                        this%sol_dlnvc(j) = (1.d0 / 1.d0 - gamma) * (this%sol_dlnvp(j) - gamma * this%sol_dlnvs(j))
                     end if
                     if (this%do_xi) then
                        this%sol_xi(j) = this%sol_vsh(j)**2/this%sol_vsv(j)**2
                     end if
                     if (this%do_ph) then
                        this%sol_ph(j) = this%sol_vpv(j)**2/this%sol_vph(j)**2
                     end if
                  end do
                  k=k+inmesh%blocks_per_layer(i)
               enddo
               
               ! convert to percental variations
               this%sol_dlnvs = this%sol_dlnvs * 100.d0
               this%sol_dlnvp = this%sol_dlnvp * 100.d0
               this%sol_dlnvsh = this%sol_dlnvsh * 100.d0
               this%sol_dlnvsv = this%sol_dlnvsv * 100.d0
               this%sol_dlnvph = this%sol_dlnvph * 100.d0
               this%sol_dlnvpv = this%sol_dlnvpv * 100.d0
               
            endif             
          
          end subroutine reparam_solution
    !========================================================



    !========================================================
    !
    !  Helper function, facilitating model export
    !
          subroutine export_solution(this,inopts,inmatr,inmesh,insche,irun,idsyn)

            implicit none
            class(post) :: this          
            class(matr) :: inmatr
            class(opts) :: inopts
            class(mesh) :: inmesh
            class(sche) :: insche

            PetscInt, intent(in) :: irun
            character(len=*), intent(in), optional :: idsyn

            PetscChar(512) ident_base
            PetscChar(512) ident

            ident_base=trim(inopts%projid)//'_'//&
                 trim(int2str(irun))//'.'//&
                 trim(int2str(int(insche%dloop(1,irun))))//'.'//&
                 trim(int2str(int(insche%dloop(2,irun))))//'.'//&
                 trim(int2str(int(insche%dloop(3,irun))))//'.'//&
                 trim(int2str(int(insche%dloop(4,irun))))

            if (present(idsyn)) then
               ident_base=trim(inopts%projid)//'_'//trim(idsyn)
               call PetscPrintf(PETSC_COMM_WORLD,"    exporting synthetic input model!\n",ierr)              
            end if
            
            call PetscPrintf(PETSC_COMM_WORLD,"    exporting model!\n",ierr)              

            if ( inmatr%processor == 0 ) then
               select case (inopts%format)
               case ('xdmf')
                  if (this%do_vsh) then
                     ident=trim(ident_base)//'_dlnvs'
                     call this%dump_model_xdmf(inmesh,insche,this%sol_dlnvs,ident)
                     ident=trim(ident_base)//'_vs'
                     call this%dump_model_xdmf(inmesh,insche,this%sol_vs,ident)
                  end if
                  if (this%do_vsv) then
                     ident=trim(ident_base)//'_dlnvp'     
                     call this%dump_model_xdmf(inmesh,insche,this%sol_dlnvp,ident)
                     ident=trim(ident_base)//'_vp'
                     call this%dump_model_xdmf(inmesh,insche,this%sol_vp,ident)  
                  end if
                  if (this%do_vc) then
                     ident=trim(ident_base)//'_dlnvc'
                     call this%dump_model_xdmf(inmesh,insche,this%sol_dlnvc,ident)
                  end if
                  if (this%do_vsh) then
                     ident=trim(ident_base)//'_dlnvph'
                     call this%dump_model_xdmf(inmesh,insche,this%sol_dlnvph,ident)
                     ident=trim(ident_base)//'_vph'
                     call this%dump_model_xdmf(inmesh,insche,this%sol_vph,ident)
                  end if
                  if (this%do_vsv) then
                     ident=trim(ident_base)//'_dlnvpv'     
                     call this%dump_model_xdmf(inmesh,insche,this%sol_dlnvpv,ident)
                     ident=trim(ident_base)//'_vpv'
                     call this%dump_model_xdmf(inmesh,insche,this%sol_vpv,ident)  
                  end if
                  if (this%do_vsh) then
                     ident=trim(ident_base)//'_dlnvsh'
                     call this%dump_model_xdmf(inmesh,insche,this%sol_dlnvsh,ident)
                     ident=trim(ident_base)//'_vsh'
                     call this%dump_model_xdmf(inmesh,insche,this%sol_vsh,ident)
                  end if
                  if (this%do_vsv) then
                     ident=trim(ident_base)//'_dlnvsv'     
                     call this%dump_model_xdmf(inmesh,insche,this%sol_dlnvsv,ident)
                     ident=trim(ident_base)//'_vsv'
                     call this%dump_model_xdmf(inmesh,insche,this%sol_vsv,ident)  
                  end if
                  if (this%do_vsv) then
                     ident=trim(ident_base)//'_xi'     
                     call this%dump_model_xdmf(inmesh,insche,this%sol_xi,ident)
                     ident=trim(ident_base)//'_ph'
                     call this%dump_model_xdmf(inmesh,insche,this%sol_ph,ident)  
                  end if
               case ('ascii')
                  if (this%do_vp) then
                     ident=trim(ident_base)//'_dlnvp.pcn'                     
                     call this%dump_model_ascii(this%sol_dlnvp,ident)
                     ident=trim(ident_base)//'_vp.abs'
                     call this%dump_model_ascii(this%sol_vp,ident)
                  end if
                  if (this%do_vs) then
                     ident=trim(ident_base)//'_dlnvs.pcn'
                     call this%dump_model_ascii(this%sol_dlnvs,ident)      
                     ident=trim(ident_base)//'_vs.abs'
                     call this%dump_model_ascii(this%sol_vs,ident)
                  end if
                  if (this%do_vc) then
                     ident=trim(ident_base)//'_dlnvc.pcn'                     
                     call this%dump_model_ascii(this%sol_dlnvc,ident)
                  end if
                  if (this%do_vsh) then
                     ident=trim(ident_base)//'_dlnvsh.pcn'
                     call this%dump_model_ascii(this%sol_dlnvsh,ident)       
                     ident=trim(ident_base)//'_vsh.abs'
                     call this%dump_model_ascii(this%sol_vsh,ident)  
                  end if
                  if (this%do_vsv) then
                     ident=trim(ident_base)//'_dlnvsv.pcn'     
                     call this%dump_model_ascii(this%sol_dlnvsv,ident)                  
                     ident=trim(ident_base)//'_vsv.abs'
                     call this%dump_model_ascii(this%sol_vsv,ident)  
                  end if
                  if (this%do_vph) then
                     ident=trim(ident_base)//'_dlnvph.pcn'
                     call this%dump_model_ascii(this%sol_dlnvph,ident)         
                     ident=trim(ident_base)//'_vph.abs'
                     call this%dump_model_ascii(this%sol_vph,ident)  
                  end if
                  if (this%do_vpv) then
                     ident=trim(ident_base)//'_dlnvpv.pcn'
                     call this%dump_model_ascii(this%sol_dlnvpv,ident)                  
                     ident=trim(ident_base)//'_vpv.abs'
                     call this%dump_model_ascii(this%sol_vpv,ident)  
                  end if
                  if (this%do_xi) then
                     ident=trim(ident_base)//'_xi.abs'
                     call this%dump_model_ascii(this%sol_xi,ident)  
                  end if
                  if (this%do_ph) then
                     ident=trim(ident_base)//'_ph.abs'
                     call this%dump_model_ascii(this%sol_ph,ident)  
                  end if 
               case default
                  call PetscPrintf(PETSC_COMM_WORLD,&
                       "POSTPROC: Unknown format, quitting!\n",ierr)
                  stop
               end select

               if (inopts%gvarr) then
                  if (.not.present(idsyn)) then
                     ident=trim(ident_base)//'_varr.dat'
                     call this%dump_varr(insche,ident,irun)
                  end if
               end if
               
            end if
          
          end subroutine export_solution
    !========================================================



    !========================================================
    !   
    !  dumps solution in a simple ascii format
    !                
          subroutine dump_run_info(this,insche,inopts)

            implicit none           
            class(post) :: this 
            class(sche) :: insche
            class(opts) :: inopts

            PetscInt i
            
            open(80,file="./results/"//trim(inopts%projid)//"_run.info")
            do i=1,insche%loop_total
               write(80,"(1x,i3,i8,2e15.7,2e15.7,2e15.7,2e15.7,2e15.7)") &
                    i,this%its(i),this%rough_abs(i),this%rough_nrm(i),&
                      this%norm_rms(i),this%norm_abs(i),this%varr_cumm(i)
            end do
            close(80)
          
          end subroutine dump_run_info
    !========================================================

    !========================================================
    !   
    !  dumps solution in a simple ascii format
    !                
          subroutine dump_varr(this,insche,ident,irun)

            implicit none           
            class(post) :: this 
            class(sche) :: insche

            character(len=*), intent(in) :: ident
            PetscInt, intent(in) :: irun
            PetscInt imat,stat
            
            open(70,file="./results/"//trim(ident),iostat=stat)
            do imat=1,insche%mats_total
               write(70,"(1x,i8,i8,i8,2e15.7)") &
                    imat,insche%fromto_info(imat,1),insche%fromto_info(imat,2),&
                    this%varr(imat,irun)
            end do
            close(70)
          
          end subroutine dump_varr
    !========================================================


    !========================================================
    !   
    !  dumps solution in a simple ascii format
    !                
          subroutine dump_model_ascii(invec,ident)

            implicit none           
            PetscScalar,      intent(in) :: invec(:)
            character(len=*), intent(in) :: ident
            PetscInt npars,stat,i
            
            open(90,file="./results/"//trim(ident),iostat=stat)
            npars = size(invec,1)
            do i=1,npars               
               write(90,"(1x,i8,2e15.7)") i,invec(i)
            end do
            close(90)
          
          end subroutine dump_model_ascii
    !========================================================



    !========================================================
    !   
    !  dumps solution in xdmf format for Paraview
    !                
          subroutine dump_model_xdmf(this,inmesh,insche,invec,ident)
                       
            implicit none           

            class(post) :: this
            class(mesh) :: inmesh
            class(sche) :: insche
            
            PetscScalar,      intent(in) :: invec(:)         
            character(len=*), intent(in) :: ident

            PetscScalar rmin,rmax           
            PetscChar(256) xdmf_elem_type
            PetscInt iinput_xdmf
            PetscInt iinput_heavy_data
            PetscInt nvertices_per_block
            PetscInt iblock,ivertex

            PetscInt i,j,k

            xdmf_elem_type = 'Hexahedron'
            nvertices_per_block = 8

            !
            !    For all unstructured topologies there is a default node 
            !    ordering. For example a HEXAHEDRON is ordered like this:
            !
            !       7 --------- 6
            !      /           /|
            !     4 --------- 5 2
            !     |  3        | /
            !     | /         |/
            !     0 --------- 1
            !        


            ! Read a stored grid from disk
            ivertex = 0
            iblock  = 0

            do i = 1,inmesh%nlays               
               rmax = r_earth - insche%layer(i)
               rmin = r_earth - insche%layer(i+1)               
               do j = 1,inmesh%blocks_per_layer(i)
                  inmesh%radmin(j,i) = rmax
                  inmesh%radmax(j,i) = rmin                  
                  k = 1 ! vertex 1
                  iblock  = iblock  + 1
                  this%connectivity(k,iblock) = ivertex + k - 1
                  this%vertices(:,ivertex+k) = sph2cart(inmesh%xlamin(j,i),&
                       inmesh%xlomin(j,i)-180.d0,inmesh%radmin(j,i))
                  k = k + 1 ! vertex 2
                  this%connectivity(k,iblock) = ivertex + k - 1
                  this%vertices(:,ivertex+k) = sph2cart(inmesh%xlamin(j,i),&
                       inmesh%xlomax(j,i)-180.d0,inmesh%radmin(j,i))
                  k = k + 1 ! vertex 3
                  this%connectivity(k,iblock) = ivertex + k - 1
                  this%vertices(:,ivertex+k) = sph2cart(inmesh%xlamax(j,i),&
                       inmesh%xlomax(j,i)-180.d0,inmesh%radmin(j,i))
                  k = k + 1 ! vertex 4
                  this%connectivity(k,iblock) = ivertex + k - 1
                  this%vertices(:,ivertex+k) = sph2cart(inmesh%xlamax(j,i),&
                       inmesh%xlomin(j,i)-180.d0,inmesh%radmin(j,i))
                  k = k + 1 ! vertex 5
                  this%connectivity(k,iblock) = ivertex + k - 1
                  this%vertices(:,ivertex+k) = sph2cart(inmesh%xlamin(j,i),&
                       inmesh%xlomin(j,i)-180.d0,inmesh%radmax(j,i))
                  k = k + 1 ! vertex 6
                  this%connectivity(k,iblock) = ivertex + k - 1
                  this%vertices(:,ivertex+k) = sph2cart(inmesh%xlamin(j,i),&
                       inmesh%xlomax(j,i)-180.d0,inmesh%radmax(j,i))
                  k = k + 1 ! vertex 7
                  this%connectivity(k,iblock) = ivertex + k - 1
                  this%vertices(:,ivertex+k) = sph2cart(inmesh%xlamax(j,i),&
                       inmesh%xlomax(j,i)-180.d0,inmesh%radmax(j,i))
                  k = k + 1 ! vertex 8
                  this%connectivity(k,iblock) = ivertex + k - 1
                  this%vertices(:,ivertex+k) = sph2cart(inmesh%xlamax(j,i),&
                       inmesh%xlomin(j,i)-180.d0,inmesh%radmax(j,i))
                  ivertex = ivertex + k
               end do
            end do

            ! xml header
            open(newunit=iinput_xdmf, file="./results/"//trim(ident)//'.xdmf')

            ! start xdmf file, write header
            write(iinput_xdmf, 733) inmesh%blocks_per_param, &
                 nvertices_per_block, 'binary', &
                 trim(ident)//'_grid.dat', &
                 ivertex, 'binary', &
                 trim(ident)//'_points.dat'

            ! create new snapshot in the temporal collection
            write(iinput_xdmf, 7341) 'grid', dble(1.d0), &
                 trim(xdmf_elem_type), inmesh%blocks_per_param, &
                 "'", "'", "'", "'"               
            
            ! write attribute
            write(iinput_xdmf, 7342) trim(ident), &
                 inmesh%blocks_per_param, 0, inmesh%blocks_per_param, 1, &
                 inmesh%blocks_per_param, trim(ident)//'_data.dat'
            write(iinput_xdmf, 7343)
            
            ! finish xdmf file
            write(iinput_xdmf, 736)
            close(iinput_xdmf)


733 format(&    
    '<?xml version="1.0" ?>',/&
    '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>',/&
    '<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">',/&
    '<Domain>',/,/&
    '<DataItem Name="grid" Dimensions="',i10, i3, '" NumberType="Int" Format="',A,'">',/&
    '  ', A,/&
    '</DataItem>',/,/&
    '<DataItem Name="points" Dimensions="',i10,' 3" NumberType="Float" Format="',A,'">',/&
    '  ', A,/&
    '</DataItem>',/,/&
    '<Grid Name="CellsTime" GridType="Collection" CollectionType="Temporal">',/)

7341 format(&    
    '    <Grid Name="', A,'" GridType="Uniform">',/&
    '        <Time Value="',F8.2,'" />',/&
    '        <Topology TopologyType="', A, '" NumberOfElements="',i10,'">',/&
    '            <DataItem Reference="/Xdmf/Domain/DataItem[@Name=', A,'grid', A,']" />',/&
    '        </Topology>',/&
    '        <Geometry GeometryType="XYZ">',/&
    '            <DataItem Reference="/Xdmf/Domain/DataItem[@Name=', A,'points', A,']" />',/&
    '        </Geometry>')

7342 format(&    
    '        <Attribute Name="', A,'" AttributeType="Scalar" Center="Cell">',/&
    '            <DataItem ItemType="HyperSlab" Dimensions="',i10,'" Type="HyperSlab">',/&
    '                <DataItem Dimensions="3 2" Format="XML">',/&
    '                    ', i10,'          0 ',/&
    '                             1          1 ',/&
    '                             1 ', i10,/&
    '                </DataItem>',/&
    '                <DataItem Dimensions="', i10, i10, '" NumberType="Float" Format="binary">',/&
    '                   ', A,/&
    '                </DataItem>',/&
    '            </DataItem>',/&
    '        </Attribute>')

7343 format(&    
    '    </Grid>',/)

736 format(&    
    '</Grid>',/&
    '</Domain>',/&
    '</Xdmf>')


            ! dump vertex data
            open(newunit=iinput_heavy_data, file="./results/"//trim(ident)//'_points.dat', access='stream', &
                 status='replace', form='unformatted', convert='little_endian')
            write(iinput_heavy_data) real(this%vertices, kind=sp)
            close(iinput_heavy_data)
            
            ! dump connectivity data
            open(newunit=iinput_heavy_data, file="./results/"//trim(ident)//'_grid.dat', access='stream', &
                 status='replace', form='unformatted', convert='little_endian')
            write(iinput_heavy_data) this%connectivity
            close(iinput_heavy_data)
            
            ! dump cell data
            open(newunit=iinput_heavy_data, file="./results/"//trim(ident)//'_data.dat', access='stream', &
                 status='replace', form='unformatted', convert='little_endian')
            write(iinput_heavy_data) real(invec, kind=sp)           
            close(iinput_heavy_data)

          
          end subroutine dump_model_xdmf
    !========================================================



    !========================================================
          function sph2cart(lat,lon,rad) result(vertex)
            
            implicit none
            PetscScalar, intent(in) :: lat,lon,rad

            PetscScalar vertex(3) ! vertex given as x,y,z
            PetscScalar phi,theta

            vertex(1) = rad * dcos(deg2rad*lat) * dcos(deg2rad*lon)
            vertex(2) = rad * dcos(deg2rad*lat) * dsin(deg2rad*lon)
            vertex(3) = rad * dsin(deg2rad*lat)

          end function sph2cart
    !========================================================                       



    !========================================================
          subroutine compute_norm(this,inmatr,irun)

            implicit none
            class(post) :: this
            class(matr) :: inmatr
            PetscInt irun

            ! Compute RMS and absolute norm
            call PetscPrintf(PETSC_COMM_WORLD,"    computing model norm!\n",ierr)
            call VecNorm(inmatr%x,NORM_2,this%norm_abs(irun),ierr)            
            call VecNorm(inmatr%x,NORM_1,this%norm_rms(irun),ierr)

            if ( verbosity > 1 ) print*,"Norm: ",this%norm_rms(irun)
            
          end subroutine compute_norm
    !========================================================


    !========================================================
          subroutine save_iterations(this,insolv,irun)

            implicit none
            class(post) :: this
            class(solv) :: insolv            
            PetscInt irun

            call PetscPrintf(PETSC_COMM_WORLD,"    storing # of iterations!\n",ierr)
            call KSPGetIterationNumber(insolv%ksp,this%its(irun),ierr)
            if ( verbosity > 1 ) print*,"Iterations: ",this%its(irun)

          end subroutine save_iterations
    !========================================================


    !========================================================
          subroutine compute_rough(this,inmatr,irun)

            implicit none
            class(post) :: this
            class(matr) :: inmatr     

            PetscInt irun

            call PetscPrintf(PETSC_COMM_WORLD,"    computing model roughness!\n",ierr)

            ! Compute adotx
            call MatMult(inmatr%A,inmatr%x,inmatr%b_rough,ierr)

            ! Apply mask vector mask_dat 
            call VecPointwiseMult(inmatr%b_rough,inmatr%b_rough,inmatr%mask_dat,ierr)

            ! Compute roughness
            call VecDot(inmatr%b_rough,inmatr%b_rough,this%rough_abs(irun),ierr)
            this%rough_nrm(irun)=this%rough_abs(irun)/&
                 this%norm_abs(irun)

          end subroutine compute_rough
    !========================================================


    !========================================================
          subroutine compute_varr(this,inmatr,insche,row_from,row_to,varr_cumm)

            implicit none
            class(post) :: this
            class(matr) :: inmatr
            class(sche) :: insche

            PetscInt, intent(in) :: row_from
            PetscInt, intent(in) :: row_to

            PetscScalar, intent(out) :: varr_cumm

            PetscScalar sum_b_sqr
            PetscScalar sum_b_min_b_adotx_sqr

            PetscInt j
            
            ! Define mask vector
            if ( inmatr%processor == 0 ) then
               do j=1,insche%rows_total+ &
                      insche%nrdamprows
                  ! Initialize mask vector with zeros
                  call VecSetValue(inmatr%mask_sub,j-1,0.d0,INSERT_VALUES,ierr)
               end do
               do j=row_from,row_to
                  ! Set mask vector 1.0 where submatrice rows are
                  call VecSetValue(inmatr%mask_sub,j,1.d0,INSERT_VALUES,ierr)
               end do
            end if

            call inmatr%assemble_vectors()
            call MatMult(inmatr%A,inmatr%x,inmatr%b_adotx,ierr)
            call VecPointwiseMult(inmatr%b_dummy,inmatr%b,inmatr%mask_sub,ierr)            
            call VecPointwiseMult(inmatr%b_adotx,inmatr%b_adotx,inmatr%mask_sub,ierr)            
            call VecDot(inmatr%b_dummy,inmatr%b_dummy,sum_b_sqr,ierr)
            call VecAXPY(inmatr%b_adotx,-1.d0,inmatr%b_dummy,ierr) 

            call VecDot(inmatr%b_adotx,inmatr%b_adotx,sum_b_min_b_adotx_sqr,ierr)
            varr_cumm = 1.d0 - (sum_b_min_b_adotx_sqr/sum_b_sqr)
            
          end subroutine compute_varr
    !========================================================



      end module module_postproc
!============================================================






!============================================================
! 
!  main program, calls all submodules
!
   program main

     use global_param

     use module_postproc
     use module_schedule
     use module_options
     use module_matrix
     use module_solver
     use module_bgmod
     use module_mesh

     implicit none

     PetscScalar outdmy

     type(opts) my_opts ! read the inversion options
     type(mesh) my_mesh ! setup the mesh
     type(sche) my_sche ! read the inversion schedule
     type(matr) my_matr ! setup the global system
     type(solv) my_solv ! setup solver context
     type(post) my_post ! background model type

      
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!           Beginning of main program routine program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      ! Initialize PETsc -> may go to a subroutine
      call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
      call PetscMemorySetGetMaximumUsage(ierr) 

      ! Welcome message
      call PetscPrintf(PETSC_COMM_WORLD,"\n                      _|                          _|\n",ierr)
      call PetscPrintf(PETSC_COMM_WORLD,"_|_|_|      _|_|    _|_|_|_|    _|_|_|    _|_|_|      _|_|_|    _|      _|\n",ierr)
      call PetscPrintf(PETSC_COMM_WORLD,"_|    _|  _|_|_|_|    _|      _|_|      _|        _|  _|    _|  _|      _|\n",ierr)
      call PetscPrintf(PETSC_COMM_WORLD,"_|    _|  _|          _|          _|_|  _|        _|  _|    _|    _|  _|\n",ierr) 
      call PetscPrintf(PETSC_COMM_WORLD,"_|_|_|      _|_|_|      _|_|  _|_|_|      _|_|_|  _|  _|    _|      _|\n",ierr)
      call PetscPrintf(PETSC_COMM_WORLD,"_|\n",ierr)                                                                          
      call PetscPrintf(PETSC_COMM_WORLD,"_|\n\n",ierr)
      call PetscPrintf(PETSC_COMM_WORLD,"Author: Ludwig Auer (ludwig.auer@gmail.com)\n\n",ierr)

      ! Read command line arguments
      call my_opts % parse_options()  
      call my_opts % tokenize_param()

      ! Initialize my mesh on basis of the command line arguments
      call PetscPrintf(PETSC_COMM_WORLD,"\n==============================\n",ierr)
      call PetscPrintf(PETSC_COMM_WORLD,"STARTING MESHER \n",ierr)
      call PetscPrintf(PETSC_COMM_WORLD,"==============================\n\n",ierr)
      call my_mesh % setup_mesh(my_opts)

      ! Initialize schedule of matrices to read and damping scheme
      call PetscPrintf(PETSC_COMM_WORLD,"\n==============================\n",ierr)
      call PetscPrintf(PETSC_COMM_WORLD,"STARTING SCHEDULER \n",ierr)
      call PetscPrintf(PETSC_COMM_WORLD,"==============================\n\n",ierr)
      call my_sche % initialize_inversion(my_opts,my_mesh)     
      
      ! When setting up the matrix we need to know all inversion
      ! options, the schedule of submatrices to be read in and
      ! the parameters associated with my mesh, first initialize
      ! the global rectangular kernel matrix
      call PetscPrintf(PETSC_COMM_WORLD,"\n==============================\n",ierr)
      call PetscPrintf(PETSC_COMM_WORLD,"SETTING UP MATRIX\n",ierr)
      call PetscPrintf(PETSC_COMM_WORLD,"==============================\n\n",ierr)
      call my_matr % initialize_matrix(my_sche,my_mesh) 

      ! read all submatrices from disk
      call my_matr % read_submatrices(my_sche)       

      !********************************
      ! Major loop over damping values
      do irun=1,my_sche%loop_total 

         ! Display solver run:
         call PetscPrintf(PETSC_COMM_WORLD,"\n==============================\n",ierr)
         call PetscPrintf(PETSC_COMM_WORLD,"RUNNING INVERSION # "&
              //trim(int2str(irun))//" of "//trim(int2str(my_sche%loop_total))//"\n",ierr)
         call PetscPrintf(PETSC_COMM_WORLD,"of type "//trim(my_opts%type)//"\n",ierr)
         call PetscPrintf(PETSC_COMM_WORLD,"==============================\n\n",ierr)

         ! In case this is not the first irun, restore rhs
         if (irun>1) call my_matr % restore_rhs()

         ! Apply damping and parameter-scaling if requested
         call PetscPrintf(PETSC_COMM_WORLD,"\n--- APPLYING DAMPING + SCALING --->\n",ierr)

         ! @TODO: Defer this to a subroutine, ugly to do it manually
         my_matr%row = my_sche%rows_total ! set running index back to nr of data

         ! Apply damping and parameter-scaling if requested
         if (my_sche%irdamp) call my_matr % apply_rdamp(my_opts,my_mesh,my_sche,irun) 
         if (my_sche%indamp) call my_matr % apply_ndamp(my_opts,my_mesh,my_sche,irun) 
         if (my_sche%iddamp) call my_matr % apply_ddamp(my_opts,my_mesh,my_sche,irun) 
         if (my_sche%iscale) call my_matr % apply_scale(my_opts,my_mesh,my_sche,irun)          

         ! @TODO: Defer this to a subroutine, ugly to do it manually
         my_matr%row = my_sche%rows_total ! doing this twice since initializing
                                          ! postproc may modify it

         ! Assemble the matrix
         call PetscPrintf(PETSC_COMM_WORLD,"\n--- ASSEMBLING MATRIX --->\n",ierr)
         call my_matr % assemble_matrix()

         ! In case this is the first run, store rhs
         if (irun==1) call my_matr % store_rhs()   

         ! Compute synthetic rhs vector
         if (my_sche%isynth) then
            call PetscPrintf(PETSC_COMM_WORLD,"\n--- COMPUTING SYNTHETICS --->\n",ierr)
            ! In case of first run, dump synthetic input model
            if (irun==1) then               
               call my_matr % read_synth_model(my_opts,my_mesh,my_sche) ! @ TODO: why do I re-read the synthetic model
               call my_post % initialize_postproc(my_opts,my_matr,my_mesh,my_sche,my_matr%x_synth) ! @ TODO: Can't i just initialize it at one loc?
               call my_post % reparam_solution(my_opts,my_matr,my_mesh,my_sche) ! reparameterize input synthetic model
               call my_post % export_solution(my_opts,my_matr,my_mesh,my_sche,irun,'synth_input') ! dump input synthetic model
            end if
            call my_matr % assemble_matrix() ! Need to reassemble matrices and vectors
            call my_matr % compute_adotx(my_matr%x_synth,my_matr%b,my_matr%mask_dmp) ! compute rhs vector
         end if

         ! Solve the system
         call PetscPrintf(PETSC_COMM_WORLD,"\n--- SOLVING Ax=b --->\n",ierr)
         call my_solv % solve_system(my_matr,my_opts)

         ! Postprocessing: reparameterize and compute VARRED and ROUGHNESS
         call PetscPrintf(PETSC_COMM_WORLD,"\n--- POSTPROCESSING --->\n",ierr)

         ! First time reqiures recomputing of the damping matrix, and takes longer
         call my_post % initialize_postproc(my_opts,my_matr,my_mesh,my_sche)

         ! Extract some information for inversion run
         call my_post % save_iterations(my_solv,irun)
         call my_post % compute_norm(my_matr,irun)

         ! Compute model roughness
         call my_matr % apply_rdamp(my_opts,my_mesh,my_sche,0)  ! irun 0 to get unweighted operator
         call my_matr % assemble_matrix() ! Need to reassemble matrices and vectors
         call my_post % compute_rough(my_matr,irun)

         ! Compute global variamce reduction
         call my_post % compute_varr(my_matr,my_sche,0,&
              my_sche%rows_total,my_post%varr_cumm(irun))

         ! Compute grouped variance reductions
         if (my_opts%gvarr) then
            call PetscPrintf(PETSC_COMM_WORLD,"    computing grouped varr!\n",ierr)            
            do imat=1,my_sche%mats_total
               call my_post % compute_varr(my_matr,my_sche,&
                                          my_sche%fromto_info(imat,1),&
                                          my_sche%fromto_info(imat,2),&
                                          my_post%varr(imat,irun))
            end do
         end if
         
         ! Reparameterize and export solution
         call my_post % reparam_solution(my_opts,my_matr,my_mesh,my_sche)
         call my_post % export_solution(my_opts,my_matr,my_mesh,my_sche,irun) ! also dums varr

      end do
      !*******************************

      call PetscPrintf(PETSC_COMM_WORLD,"\n==============================\n",ierr)
      call PetscPrintf(PETSC_COMM_WORLD,"FINALIZING INVERSION\n",ierr)
      call PetscPrintf(PETSC_COMM_WORLD,"==============================\n\n",ierr)
      if ( my_matr%processor == 0) call my_post % dump_run_info(my_sche,my_opts)
      
      ! Free work space.
      call my_post % destroy_postproc()         
      call my_solv % destroy_solver()
      call my_matr % destroy_matrix()

      if (verbosity > 1) then
         ! Display amount of memory used
         call PetscMemoryGetMaximumUsage(memory,ierr)
         memory=memory*1.0d-6
         write(*,*),'Maximum Memory usage = ',memory,'Mb'
      end if

      ! Finalize PETSc
      call PetscFinalize(ierr)

    end program main




