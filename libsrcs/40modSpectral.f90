module Spectral
  use Globals
  use GaussQuadrature
  use DenseMatrix
  use QuadGrid
  implicit none
  private
  
  type, public :: spct
     !> Max degree of 1d polynomial
     integer :: ndeg = 0
     !> Non zero term of stiff matrix in 
     !> up. triang. form = (((ndeg+1)**2)*((ndeg+1)**2-1))/2
     integer :: nterm = 0
     !> Dimension((ndeg+1),(ndeg+1))
     integer, allocatable :: degx(:,:)
     !> Dimension((ndeg+1),(ndeg+1))
     integer, allocatable :: degy(:,:)
     !> Dimension((ndeg+1)**2)     
     integer, allocatable :: saved_dx(:)
     !> Dimension((ndeg+1)**2)
     integer, allocatable :: saved_dy(:)
     !> Dimension((ndeg+1)**2,(ndeg+1)**2)
     integer, allocatable :: dx2d(:,:)
     !> Dimension((ndeg+1)**2,(ndeg+1)**2)
     integer, allocatable :: dy2d(:,:)
     !> Dimension( (ndeg+1)**2 )
     real(kind=double), allocatable :: vander(:)
     !> Dimension( (ndeg+1)**2 )
     real(kind=double), allocatable :: integral_base(:)
   contains
     !> static constructor
     !> (procedure public for type degrees)
     procedure, public, pass :: init => init_spec
     !> static destructor
     !> (procedure public for type degrees)
     procedure, public, pass :: kill => kill_spec
     !> Info procedure.
     !> (public procedure for type degrees)
     procedure, public, pass :: info => info_spec
     !> Procedure for evaluation  degx(i,j) and degy(i,j)
     !> of type spec.
     procedure, private, nopass :: degs
     !> Procedure for building the member degx and degy
     !> of type spec. Initiliazed the members.
     !> Used in init
     procedure, private, pass :: build_degs
     !> Procedure for evaluation  dx(k) of type spec 
     procedure, public, nopass :: dx
     !> Procedure for evaluation  dy(k) of type spec 
     procedure, public, nopass :: dy
     !> Procedure for building the member dx and dy
     !> of type spec. Initiliazed the members.
     !> Used in init
     procedure, private, pass :: build_dxdy
     !> Procedure for building the member der1d, dx2d
     !> and 2y2d of type spec. Initiliazed the members.
     !> Used in init
     procedure, private, pass :: build_der 
     !> Procedure to calculate matrix CTP, CUP 
     !> for interval [a,b]
     procedure, public, pass :: eval_CTPCUP2 
     !> Procedure evalute dx_coeff, given coeff 
     !> for interval [a,b]
     procedure, public, pass :: make_dx
     !> Procedure evalute dy_coeff, given coeff
     !> for interval [a,b]
     procedure, public, pass :: make_dy
     !> Procedure to compute function on point (x,y) 
     procedure, public, pass :: eval
     !> Procedure x component of the grad on point (x,y) 
     procedure, public, pass :: grad_x
     !> Procedure y component of the grad on point (x,y)
     procedure, public, pass :: grad_y
     !> Procedure nrm of gradient on a point(x,y) 
     procedure, public, pass :: nrm_grad
     !> Procedure evalute dx_coeff, given coeff 
     !> for interval [a,b]
     procedure, public, nopass :: eval_vander1d
     !> Procedure evalute dx_coeff, given coeff 
     !> for interval [a,b]
     procedure, public, pass :: make_vander1d
     !> Procedure evalute dy_coeff, given coeff
     !> for interval [a,b]
     procedure, public, pass :: make_vander1d_points
     !> Procedure evalute dy_coeff, given coeff
     !> for interval [a,b]
     procedure, public, pass :: vander2d
     !> Procedure evalute dy_coeff, given coeff
     !> for interval [a,b]
     procedure, public, pass :: make_vander2d_points
     !> Procedure 
     procedure, public, pass :: make_vander_grid
     !> Procedure 
     procedure, public, pass :: make_vander
  end type spct

  type, public :: inter
     !> Interval extremes
     real(kind=double) :: a, b
     !> Dimension (spec%ndeg+1,gauss_inter%ngauss)
     real(kind=double), allocatable :: vander_gauss(:,:)
     !> Dimension ( gauss_inter%ngauss )
     real(kind=double), allocatable :: coord(:)
     !> Dimension ( gauss_inter%ngauss )
     real(kind=double), allocatable :: weight(:)
   contains
     !> static constructor
     !> (procedure public for type interval)
     procedure, public, pass :: init => init_inter
     !> static destructor
     !> (procedure public for type interval)
     procedure, public, pass :: kill => kill_inter
     !> Info procedure.
     !> (public procedure for type interval)
     procedure, public, pass :: info => info_inter
     !> Procedure for building member vander_gauss
     !> (public procedure for type interval, used in init)
     procedure, public, pass :: make_vander_gauss
  end type inter

  
contains
  !>-------------------------------------------
  !> Initizialization procedure for type spec
  !<-------------------------------------------
  subroutine init_spec(this,lun_err,ndeg)
    use Globals
    implicit none
    class(spct), intent(inout) :: this        
    integer, intent(in) :: lun_err
    integer, intent(in) :: ndeg
    !local
    logical :: rc
    integer :: res
    integer :: i,j,m

    this%ndeg  = ndeg
    this%nterm = (((ndeg+1)**2)*((ndeg+1)**2-1))/2

    call this%build_degs(lun_err)
    call this%build_dxdy(lun_err)
    call this%build_der(lun_err)
    allocate(&
         this%vander  ( (ndeg+1)**2 ),&
         this%integral_base( (ndeg+1)**2 ),&
         stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'init_spec', &
         ' type spec member coeff, dx_coeff, dy_coeff ',res)

    do m = 1,(ndeg+1)**2 
       i = dx(m,ndeg)
       j = dy(m,ndeg)
       if (( i .eq. 1) .or. (j .eq. 1) ) then
          this%integral_base(m)=zero
       else
          this%integral_base(m) = &
            one / (one-i**2) * (one+(-1)**i) * &
            one / (one-j**2) * (one+(-1)**j)
       end if
    end do
    

  end subroutine init_spec


  !>-------------------------------------------
  !> Desctruction procedure for type spec
  !<-------------------------------------------
  subroutine kill_spec(this,lun_err)
    use Globals
    implicit none
    class(spct), intent(inout) :: this        
    integer, intent(in) :: lun_err
    !local 
    logical rc
    integer res

    deallocate(&
         this%degx,&
         this%degy,&
         this%saved_dx,&
         this%saved_dy,&
         this%dx2d,&
         this%dy2d,&
         this%vander,&
         this%integral_base,&
         stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, 'main', &
         ' type spec member degx, degy, dx ,dy, dx2d, dy2d ',res)

  end subroutine kill_spec
  !>-------------------------------------------
  !> Info procedure for type spec
  !<-------------------------------------------
  subroutine info_spec(this,lun_out)
    use Globals
    implicit none
    class(spct), intent(inout) :: this        
    integer, intent(in) :: lun_out
    !local 
    logical rc
    integer res

    write(lun_out,*) 'Maximum degree of 1d polynomial=',this%ndeg

  end subroutine info_spec

  !>--------------------------------------------
  !> public function for calculatation of 
  !> degx(i,j) and degy(i,j)
  !<-------------------------------------------
  function degs(axe_label,i,j)
    implicit none
    character(len=1), intent(in) :: axe_label
    integer,          intent(in) :: i,j
    integer                      :: degs
    !local 
    if (axe_label == 'x') then
       degs = j-1
    end if
    if (axe_label == 'y') then
       degs = i-1
    end if
  end function degs

  !>------------------------------------------------
  !> Procedure for building the member degx and degy
  !> of type spec. Initiliazed the members.
  !<------------------------------------------------
  subroutine build_degs(this,lun_err)
    use Globals
    implicit none
    class(spct), intent(inout) :: this
    integer,     intent(in   ) :: lun_err
    !local 
    logical rc
    integer i,j,res
    !integer degs

    allocate(&
         this%degx(this%ndeg+1,this%ndeg+1),&
         this%degy(this%ndeg+1,this%ndeg+1),&
         stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'build_id_inter', &
         ' type mesh member id_inter',res)

    do j = 1, this%ndeg+1           
       do i = 1, this%ndeg+1
          this%degx(i,j) = degs('x',i,j)
          this%degy(i,j) = degs('y',i,j)
       end do
    end do
  end subroutine build_degs

  !>----------------------------------
  !> Procedure for the evalution of dx
  !>---------------------------------
  function dx(m,ndeg)
    use Globals
    implicit none
    integer m,ndeg
    integer dx
    dx = (m-1)/(ndeg+1)
  end function dx

  !>----------------------------------
  !> Procedure for the evalution of dy
  !>---------------------------------
  function dy(m,ndeg)
    use Globals
    implicit none
    integer m,ndeg
    integer dy
    dy = mod(m-1,ndeg+1)
  end function dy


  !>------------------------------------------------
  !> Procedure for building the member dx and dy
  !> of type spec. Initiliazed the members.
  !<------------------------------------------------
  subroutine build_dxdy(this,lun_err)
    use Globals
    implicit none
    class(spct), intent(inout) :: this
    integer,     intent(in   ) :: lun_err
    !local 
    logical rc
    integer i,j,m,res

    allocate(&
         this%saved_dx( (this%ndeg+1)**2 ),&
         this%saved_dy( (this%ndeg+1)**2 ),&
         stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'build_dxdy', &
         ' type mesh member dx, dy ',res)

    do m = 1, (this%ndeg+1)**2 
       this%saved_dx(m) = dx(m,this%ndeg)
       this%saved_dy(m) = dy(m,this%ndeg)
    end do
  end subroutine build_dxdy

  !>------------------------------------------------
  !> Procedure for building the member dx2d 
  !> and dy2d of type spec. Initiliazed the members.
  !<------------------------------------------------
  subroutine build_der(this,lun_err)
    use Globals
    implicit none
    class(spct), intent(inout) :: this
    integer,     intent(in   ) :: lun_err
    !local 
    logical rc
    integer i,j,m,res,start
    integer, allocatable :: der1d(:,:)
    integer, allocatable :: b(:)      

    allocate(&
         this%dx2d( (this%ndeg+1)**2 , (this%ndeg+1)**2 ),&
         this%dy2d( (this%ndeg+1)**2 , (this%ndeg+1)**2 ),&
         stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'build_der', &
         ' type mesh member der1d, dx2d, dx2d ',res)

    allocate(&
         der1d(this%ndeg+1,this%ndeg+1),&
         b( (this%ndeg+1)**2  ),&
         stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'build_der', &
         ' local array b ',res)

    ! make der1d
    der1d = 0
    do j = 2, this%ndeg+1
       start = mod(j,2)+1
       do i = start, j-1, 2
          der1d(i,j) = 2*(j-1)
       end do
    end do
    der1d(1,:) = der1d(1,:)/2


    ! make der2d
    this%dy2d=0
    this%dx2d=0
    do i=1,this%ndeg+1
       this%dy2d(&
            1+(i-1)*(this%ndeg+1):i*(this%ndeg+1),   &
            1+(i-1)*(this%ndeg+1):i*(this%ndeg+1)) = &
            der1d
    end do
    m=0
    do j = 1, this%ndeg+1
       do i = 1, this%ndeg+1
          m = m + 1
          b(m) = 1 + (j-1) + (this%ndeg+1)*(i-1)
       end do
    end do
    this%dx2d = this%dy2d(b,b)

    deallocate(&
         der1d,&
         b,&
         stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, 'build_der', &
         ' var local b',res)
  end subroutine build_der

!!$  !>--------------------------------------------
!!$  !> Procedure For evaluation of matrix CTP, CUP
!!$  !> on an interval [a,b]
!!$  !>--------------------------------------------
!!$  subroutine eval_CTPCUP(this,ndeg,a,b,CTP,CUP)
!!$    use Globals
!!$    implicit none
!!$    class(spct),       intent(in ) :: this
!!$    integer,           intent(in ) :: ndeg
!!$    real(kind=double), intent(in ) :: a,b
!!$    real(kind=double), intent(out) :: CTP(ndeg+1,ndeg+1)
!!$    real(kind=double), intent(out) :: CUP(ndeg+1,ndeg+1)
!!$    !local
!!$    integer :: iinter,i,j
!!$    integer :: n1,n2,n3,n4
!!$    real(kind=double) :: low_ang, up_ang
!!$
!!$    low_ang = acos(a)
!!$    up_ang  = acos(b)
!!$
!!$    do j=1,ndeg+1
!!$       do i=1,ndeg+1
!!$          n1 = this%degx(i,j) + this%degy(i,j) + 1
!!$          n2 = this%degx(i,j) + this%degy(i,j) - 1
!!$          n3 = this%degx(i,j) - this%degy(i,j) + 1
!!$          n4 = this%degx(i,j) - this%degy(i,j) - 1
!!$          CTP(i,j) = &
!!$               sint(n1,low_ang,up_ang) &
!!$               - sint(n2,low_ang,up_ang) &
!!$               + sint(n3,low_ang,up_ang) &
!!$               - sint(n4,low_ang,up_ang)
!!$       end do
!!$    end do
!!$    call assCU3(ndeg+1,ndeg+1,CTP,CUP)
!!$  contains
!!$
!!$    function sint(ngrade,low_ang,up_ang) result (si)
!!$      use Globals
!!$      integer :: ngrade
!!$      real(kind=double) :: low_ang,up_ang
!!$      real(kind=double) :: si
!!$      if (ngrade == 0) then
!!$         si=zero
!!$      else
!!$         si = - &
!!$              sin(ngrade*(up_ang+low_ang)/2.0d0) * &
!!$              sin(ngrade*(up_ang-low_ang)/2.0d0) / (2.0d0*ngrade)
!!$      end if
!!$    end function sint
!!$
!!$    
!!$  end subroutine Eval_CTPCUP
!!$
!!$  subroutine assCU3(nrow,ncol,mat_in,mat_out)
!!$      use Globals
!!$      implicit none 
!!$      integer :: nrow,ncol
!!$      real(kind=double), intent(in   ) :: mat_in(nrow,ncol)
!!$      real(kind=double), intent(out  ) :: mat_out(nrow,ncol)
!!$      ! local
!!$      integer :: i,j
!!$
!!$
!!$      mat_out = zero     
!!$      do j = 1 , ncol - 1
!!$         ! copy the column j swifted one down and one left
!!$         call dcopy(nrow-1,mat_in(1,j),1,mat_out(2,j+1),1)
!!$         ! scal the first column by onehalf
!!$         if (j .eq. 1) then
!!$            call dscal(nrow-1,onehalf,mat_out(2,j+1),1)
!!$         end if
!!$         ! scal the first row by onehalf
!!$         mat_out(2,j+1) = mat_out(2,j+1) * onehalf
!!$         do i = 3 , nrow- 1
!!$            mat_out(i+1,j+1) = mat_out(i+1-2,j+1) + mat_in(i,j)
!!$         end do
!!$         if( j .gt. 1) then
!!$            !mat_out(2:nrow,j+1) = mat_out(2:nrow,j+1) + mat_out(2:nrow,j+1-2)
!!$            call daxpy(nrow-1,1.0d0,mat_out(2,j+1-2),1,mat_out(2,j+1),1)
!!$         end if
!!$      end do
!!$      !mat_out=4.0d0*mat_out
!!$      ! mat_ou(2,2) first non zero element
!!$      call dscal((nrow-1)*ncol-1,4.0d0,mat_out(2,2),1)
!!$
!!$    end subroutine assCU3

  !>--------------------------------------------
  !> Procedure For evaluation of matrix CTP, CUP
  !> on an interval [a,b]
  !>--------------------------------------------
  subroutine eval_CTPCUP2(this,ndeg,a,b,CTP,CUP)
    use Globals
    implicit none
    class(spct),       intent(in ) :: this
    integer,           intent(in ) :: ndeg
    real(kind=double), intent(in ) :: a,b
    real(kind=double), intent(out) :: CTP(ndeg+1,ndeg+1)
    real(kind=double), intent(out) :: CUP(ndeg+1,ndeg+1)
    !local
    integer :: iinter,i,j
    integer :: n1,n2,n3,n4
    real(kind=double) :: low_ang, up_ang, sum , dif

    if ( ( a > b) .or. ( abs(a) > one ) .or. ( abs(b)> one ) ) then
       write(*,*) 'Error in eval_CTPCUP2, not andmissible intervals'
    end if
    low_ang = acos(a)
    up_ang  = acos(b)
    sum = up_ang + low_ang
    dif = up_ang - low_ang
    do j=1,ndeg+1
       do i=1,ndeg+1
          n1 = this%degx(i,j) + this%degy(i,j) + 1
          n2 = this%degx(i,j) + this%degy(i,j) - 1
          n3 = this%degx(i,j) - this%degy(i,j) + 1
          n4 = this%degx(i,j) - this%degy(i,j) - 1
          !write(*,*) i,j,n1,n2,n3,n4
          CTP(i,j) = &
               sint(n1,sum,dif) &
               - sint(n2,sum,dif) &
               + sint(n3,sum,dif) &
               - sint(n4,sum,dif)
       end do
    end do



    call assCU3(ndeg+1,ndeg+1,CTP,CUP)

  contains

    function sint(ngrade,sum,diff) result (si)
      use Globals
      integer :: ngrade
      real(kind=double) :: sum, diff
      real(kind=double) :: si
      if (ngrade == 0) then
         si=zero
      else
         si = - &
              sin(ngrade*(sum)/2.0d0) * &
              sin(ngrade*(dif)/2.0d0) / (2.0d0*ngrade)
      end if
    end function sint

    
    
  end subroutine Eval_CTPCUP2

subroutine assCU3(nrow,ncol,mat_in,mat_out)
      use Globals
      implicit none 
      integer :: nrow,ncol
      real(kind=double), intent(in   ) :: mat_in(nrow,ncol)
      real(kind=double), intent(out  ) :: mat_out(nrow,ncol)
      ! local
      integer :: i,j
      mat_out = zero     
      do j = 1 , ncol - 1
         ! copy the column j swifted one down and one left
         call dcopy(nrow-1,mat_in(1,j),1,mat_out(2,j+1),1)

         ! scal the first column by onehalf
         if (j .eq. 1) then
            call dscal(nrow-1,onehalf,mat_out(2,j+1),1)
         end if
         ! scal the first element of the column by onehalf
         mat_out(2,j+1) = mat_out(2,j+1) * onehalf
         if  (j == 1) then
            do i = 2 , nrow - 1
               mat_out(i+1,j+1) = mat_out(i+1-2,j+1) + mat_in(i,j) * onehalf
            end do
         else
            !mat_out(i+1,j+1) = mat_out(i+1-2,j+1) + mat_in(i,j) 
            do i = 2 , nrow - 1
               mat_out(i+1,j+1) = mat_out(i+1-2,j+1) + mat_in(i,j) 
            end do
         end if
         if( j .gt. 1) then
            !mat_out(2:nrow,j+1) = mat_out(2:nrow,j+1) + mat_out(2:nrow,j+1-2)
            call daxpy(nrow-1,1.0d0,mat_out(2,j+1-2),1,mat_out(2,j+1),1)
         end if
      end do
      !mat_out=4.0d0*mat_out
      ! mat_ou(2,2) first non zero element
      call dscal((nrow-1)*ncol-1,4.0d0,mat_out(2,2),1)

    end subroutine assCU3

    !> This procedure compute the coefficient of 
    !> gradx_pot=\sum_{i} coeff_i \Psi_i(x,y)
    !> given pot=\sum_{i} dx_coeff_i \Psi (X,y)
    !> This work as Block matrix with  such that
    !>    dx_coeff = A coeff (ndeg+1)x(ndeg+1) blocks and
    !> with A=(
    !>     0     1*Id  0     3*Id  0     5*Id 0     7*Id  0    ...
    !>     0     0     4*Id  0     8*Id  0    12*Id 0    16*Id ...
    !>     0     0     0     6*Id  0    10*Id 0    14*Id  0    ...
    !>     0     0     0     0     8*Id     0 12*Id 0    16*Id ...
    !>     0     0     0     0     0    10*Id 0    14*Id  0    ...
    !>     0     0     0     0     0     0    12*Id 0    16*Id ...
    !>     0     0     0     0     0     0     0    14*Id  0   ...
    !>     0     0     0     0     0     0     0     0    16*Id...
    !>     0     0     0     0     0     0     0     0     0  ...
    !>) where Id is the Identity matrix of size (ndeg+1)x(ndeg+1)
    !TODO definitely not optimal and not clear
   subroutine make_dx(this,coeff,dx_coeff)
    use Globals
    implicit none
    class(spct),       intent(in ) :: this
    real(kind=double), intent(in ) :: coeff( (this%ndeg+1)**2 )
    real(kind=double), intent(inout) :: dx_coeff( (this%ndeg+1)**2 )
    ! local
    integer i,j, iblock, jblock, j_loc, k, ndeg

    ndeg=this%ndeg

    do i=1, (ndeg+1)**2
       dx_coeff(i) = zero 
       iblock = (i-1) / (ndeg+1) +1
       j_loc  = mod( (i-1) , ndeg+1) +1
       if (iblock == 1 ) then
          !write(*,*) iblock, jloc, 
          do k = iblock, ndeg, 2
             j = k * ( ndeg + 1) + j_loc
             dx_coeff(i) = dx_coeff(i) + (2*(k/2) + 1 ) * coeff(j)
             !write(*,*) iblock, j_loc, k,j, coeff(j), dx_coeff(i)

          end do
       else
          do k = iblock , ndeg, 2
             j = k * ( ndeg + 1) + j_loc
             dx_coeff(i) = dx_coeff(i) + 2 * (k) * coeff(j)
             !write(*,*) iblock, j_loc, k,j, coeff(j), dx_coeff(i)

          end do
       end if
    end do
  end subroutine make_dx

  !> This procedure compute the coefficient of 
  !> grady_pot=\sum_{i} coeff_i \Psi_i(x,y)
  !> given pot=\sum_{i} dy_coeff_i \Psi (X,y)
  !> This work as a block diagonal matrix A with
  !> (ndeg+1)x(ndeg+1) blocks acting on coeff, thus
  !>
  !>    dy_coeff = A coeff 
  !>
  !> A=diag(B)
  !> with B=(
  !>      0     1     0     3     0     5     0     7     0
  !>      0     0     4     0     8     0    12     0    16
  !>      0     0     0     6     0    10     0    14     0
  !>      0     0     0     0     8     0    12     0    16
  !>      0     0     0     0     0    10     0    14     0
  !>      0     0     0     0     0     0    12     0    16
  !>      0     0     0     0     0     0     0    14     0
  !>      0     0     0     0     0     0     0     0    16
  !>      0     0     0     0     0     0     0     0     0
  !>)
  !TODO definitely not optimal and not clear
  subroutine make_dy(this,coeff,dy_coeff)
    use Globals
    implicit none
    class(spct), intent(inout) :: this
    real(kind=double), intent(in ) :: coeff( (this%ndeg+1)**2 )
    real(kind=double), intent(inout) :: dy_coeff( (this%ndeg+1)**2 )

    ! local
    integer i, begin, finish

    do i = 1, this%ndeg+1
       begin  = (i-1) * (this%ndeg+1) +1
       finish = i     * (this%ndeg+1)
       call block_product(this%ndeg, &
            coeff(begin:finish), &
            dy_coeff(begin:finish))
    end do

  contains
    subroutine block_product(ndeg, coeff_block, dy_block)
      implicit none
      integer ndeg
      real(kind=double), intent(in ) :: coeff_block(ndeg+1)
      real(kind=double), intent(out) :: dy_block(ndeg+1)
      !Local
      integer k,l

      dy_block=zero
      ! 1
      do k = 2, ndeg+1, 2  
         dy_block(1) = dy_block(1) + (2*(k/2) - 1 ) * coeff_block(k)
      end do
      ! 1:ndeg
      do l = 2, ndeg
         do k = l+1, ndeg+1, 2
            dy_block(l) = dy_block(l) + 2 * (k-1) * coeff_block(k)
         end do
      end do
      ! ndeg+1
      !dy_block(ndeg+1)=zero
    end subroutine block_product

  end subroutine make_dy

  function eval_vander1d(ideg, x)
    use Globals
    integer,           intent(in) :: ideg
    real(kind=double), intent(in) :: x
    real(kind=double) :: eval_vander1d

    eval_vander1d = cos( ideg * acos(x) )

  end function eval_vander1d

  subroutine make_vander1d(this, x, vander1d)
    use Globals
    class(spct),       intent(in ) :: this
    real(kind=double), intent(in ) :: x
    real(kind=double), intent(out) :: vander1d(this%ndeg+1)
    !local
    integer :: i
    real(kind=double) :: ang
    ang = acos(x)
    do i = 1, this%ndeg+1
       vander1d(i) = cos( (i-1) * ang )
    end do
  end subroutine make_vander1d

  subroutine make_vander1d_points(this, n1dpoints, x, vander1d_points)
    use Globals
    class(spct),       intent(in ) :: this
    integer,           intent(in ) :: n1dpoints 
    real(kind=double), intent(in ) :: x(n1dpoints)
    real(kind=double), intent(out) :: vander1d_points(this%ndeg+1,n1dpoints)
    !local
    integer :: i,k
    do k=1,n1dpoints
       call this%make_vander1d( x(k), vander1d_points(:,k) )
    end do
  end subroutine make_vander1d_points

  subroutine make_vander(this,x,y)
    implicit none
    class(spct),       intent(inout) :: this
    real(kind=double), intent(in   ) :: x, y
    ! local
    integer i,j,m
    real(kind=double) :: temp
    
    !m=0
    do m = 1,(this%ndeg+1)**2 
       i = dx(m,this%ndeg)
       j = dy(m,this%ndeg)       
       this%vander(m) = this%eval_vander1d(i, x) * this%eval_vander1d(j, y)
    end do
!!$    do i = 1, this%ndeg + 1
!!$       temp = this%eval_vander1d(i, x) 
!!$       do j=1, this%ndeg  + 1
!!$          m=m+1
!!$          this%vander(m) = temp * this%eval_vander1d(j, y)
!!$       end do
!!$    end do
  end subroutine make_vander
    

  function vander2d(this, id_deg, x, y)
    use Globals
    implicit none
    class(spct),       intent(in) :: this
    integer,           intent(in) :: id_deg
    real(kind=double), intent(in) :: x, y
    !output
    real(kind=double) :: vander2d
    vander2d = &
         cos( this%dx(id_deg, this%ndeg) * acos( x ) ) * &
         cos( this%dy(id_deg, this%ndeg) * acos( y ) )
  end function vander2d

  subroutine make_vander2d_points(this, n2dpoints, coord, vander2d_points)
    use Globals
    class(spct),       intent(in ) :: this
    integer,           intent(in ) :: n2dpoints 
    real(kind=double), intent(in ) :: coord(2,n2dpoints)
    real(kind=double), intent(out) :: vander2d_points((this%ndeg+1) **2 , n2dpoints)
    !local
    integer :: i,k,m
    do k=1,n2dpoints
       do m=1, (this%ndeg +1) **2
          vander2d_points(m,k) = this%vander2d(m, coord(1,k), coord(2,k))
       end do
    end do
  end subroutine make_vander2d_points


  subroutine make_vander_grid(this,&
       nx, ny, &
       xcoord, ycoord,&
       vander_grid)
    use Globals
    class(spct),       intent(in ) :: this
    integer,           intent(in ) :: nx, ny
    real(kind=double), intent(in ) :: xcoord(nx), ycoord(ny)
    real(kind=double), intent(out) :: vander_grid((this%ndeg+1)**2 , nx*ny)
    !local
    integer :: ix,iy,k, m
    k=0
    do ix = 1, nx
       do iy = 1, ny
          k=k+1
          do m=1, (this%ndeg +1) **2
             vander_grid(m,k)= this%vander2d(m, xcoord(ix), ycoord(iy) )
          end do
       end do
    end do
  end subroutine make_vander_grid

  function eval(this,coeff,x,y)
    implicit none
    class(spct),       intent(inout) :: this
    real(kind=double), intent(in   ) :: coeff( ( this%ndeg+1 )**2 )
    real(kind=double), intent(in   ) :: x,y
    
    real(kind=double) :: eval
    real(kind=double) :: ddot
        
    call this%make_vander(x,y)
    eval = ddot( (this%ndeg+1)**2, this%vander, 1, coeff,1)
    
  end function eval

  function grad_x(this,dx_coeff,x,y)
    implicit none
    class(spct),       intent(inout) :: this
    real(kind=double), intent(in   ) :: dx_coeff( ( this%ndeg+1 )**2 )
    real(kind=double), intent(in   ) :: x,y
    
    real(kind=double) :: grad_x
    real(kind=double) :: ddot

    call this%make_vander(x,y)
    grad_x = ddot( (this%ndeg+1)**2, &
         this%vander,1,&
         dx_coeff,1)

  end function grad_x

  function grad_y(this, dy_coeff,x,y)
    implicit none
    class(spct),       intent(inout) :: this
    real(kind=double), intent(in   ) :: dy_coeff( (this%ndeg+1)**2 )
    real(kind=double), intent(in   ) :: x,y
    
    real(kind=double) :: grad_y
    real(kind=double) :: ddot

    call this%make_vander(x,y)
    grad_y = ddot((this%ndeg+1)**2,&
         this%vander,1,&
         dy_coeff,1)

  end function grad_y

  function nrm_grad(this, ndeg,dx_coeff,dy_coeff,x,y)
    implicit none
    class(spct),       intent(inout) :: this
    integer,           intent(in   ) :: ndeg

    real(kind=double), intent(in   ) :: dx_coeff( ( ndeg+1 )**2 )
    real(kind=double), intent(in   ) :: dy_coeff( ( ndeg+1 )**2 )
    real(kind=double), intent(in   ) :: x,y
    real(kind=double) :: nrm_grad
    !local
    real(kind=double) :: ddot
    real(kind=double) :: grad_x, grad_y
    call this%make_vander(x,y)
    grad_x = ( ddot( ( ndeg+1 )**2 , &
         this%vander,1,&
         dx_coeff,1) )**2 
    grad_y = ( ddot( ( ndeg+1 )**2,&
         this%vander,1,&
         dy_coeff,1) )**2 

    nrm_grad = sqrt(grad_x**2+grad_y**2)
        
  end function nrm_grad

  subroutine init_inter(this,&
       lun_err,&
       a,b, &
       spec,&
       gauss_inter)
    use Globals
    use GaussQuadrature
    implicit none
    class(inter),      intent(inout) :: this
    integer,           intent(in   ) :: lun_err
    real(kind=double), intent(in   ) :: a,b 
    type(spct),        intent(in   ) :: spec
    type(gaussq),      intent(inout) :: gauss_inter

    
    !localg
    logical :: rc
    integer :: res
    
    this%a = a
    this%b = b
    allocate (&
!!$       this%CTP( spec%ndeg+1, spec%ndeg+1 ),&
!!$       this%CUP( spec%ndeg+1, spec%ndeg+1 ),&
       this%vander_gauss(spec%ndeg+1, gauss_inter%ngauss ),&
       this%coord(gauss_inter%ngauss ),&
       this%weight(gauss_inter%ngauss ),&
       stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_alloc, 'init_inteval', &
       ' member coeff, weight',res)

!!$    call spec%eval_CTPCUP2(spec%ndeg,&
!!$               this%a,&
!!$               this%b,&
!!$               this%CTP(:,:),&
!!$               this%CUP(:,:))
    !call this%make_vander_gauss(gauss_inter, spec)
  end subroutine init_inter

  subroutine kill_inter(this,lun_err)
    use Globals
    implicit none
    class(inter), intent(inout) :: this
    integer,      intent(in   ) :: lun_err
    !local
    logical :: rc
    integer :: res
    
    this%a = zero
    this%b = zero
    deallocate (&
!!$       this%CTP,&
!!$       this%CUP,&
         this%coord,&
         this%weight,&
         this%vander_gauss,&
       stat=res)
    if(res .ne. 0) rc = IOerr(lun_err, err_dealloc, 'kill_interval', &
       ' member coord, weight, vander_gauss',res)
  end subroutine kill_inter
    

  subroutine info_inter(this,lun_out)
    use Globals
    implicit none
    class(inter), intent(in) :: this
    integer,      intent(in) :: lun_out
    !local
    integer ndeg,ngauss

    ndeg  =size(this%vander_gauss, 1)-1
    ngauss=size(this%weight)
    write(*,*) ' a=', this%a,' b=', this%b
    write(*,*)  'ndeg =',ndeg,'ngauss', ngauss
    
!!$    ndeg =  size( this%CTP,1) -1 
!!$    call print_rmat(ndeg+1,ndeg+1, this%CTP)
!!$    call print_rmat(ndeg+1,ndeg+1, this%CUP)
    


  end subroutine info_inter
  
  subroutine make_vander_gauss(this, gauss_inter, spec)
    use Globals
    use GaussQuadrature
    class(inter), intent(inout) :: this 
    type(gaussq), intent(inout) :: gauss_inter
    type(spct),   intent(in   ) :: spec
    !local
    integer :: iinter
    real(kind=double) :: a,b, len_inter

    call gauss_inter%on_interval(this%a,this%b)
    this%coord  = gauss_inter%coord_ab
    this%weight = gauss_inter%weight_ab

    call spec%make_vander1d_points(&
         gauss_inter%ngauss, gauss_inter%coord_ab, &
         this%vander_gauss(:,:) )
  end subroutine make_vander_gauss

       



end module Spectral

!!$function  integral_nrm_grad2(&
!!$     ix_inter,iy_inter,&
!!$     dx_coeff, dy_coeff,&
!!$     basis1d_ongauss_xinter,&
!!$     basis1d_ongauss_yinter,&
!!$     vander_loc,&
!!$     spec)
!!$  integer :: ndeg, ngauss, ninter
!!$  integer :: ix_inter, iy_inter
!!$  real(kind=double) :: basis1d_onguass_xinter( ndeg, ngauss)
!!$  real(kind=double) :: basis1d_onguass_yinter( ndeg, ngauss)
!!$  real(kind=double) :: vander_loc( (ndeg+1)**2 ) 
!!$  type(spct), intent(in) :: spec
!!$
!!$  real(kind=double) :: integral_nrm_grad
!!$  
!!$  integral_nrm_grad = zero
!!$  do ixguass = 1, spec%ngauss
!!$     do iygauss = 1, spec%ngauss
!!$        k = k + 1
!!$        call vander_loc(&
!!$             ndeg, ngauss, ninter&
!!$             ix_inter, iy_inter, ix_gauss, iy_gauss,&
!!$             basis1d_onguass_xinter,&
!!$             basis1d_onguass_yinter,&
!!$             vander_loc)
!!$        
!!$        integral_nrm_grad = integral_nrm_grad +  &
!!$             spec%gauss_1dweight(ixgauss) * spec%gauss1d_weigth(iy_gauss) * sqrt(&
!!$             ( ddot((ndeg+1)**2,vander_loc(1),1,dx_coeff,1) ) ** 2 + &
!!$             ( ddot((ndeg+1)**2,vander_loc(1),1,dx_coeff,1) ) ** 2 )
!!$     end do
!!$  end do
!!$  
!!$  contains
!!$    subroutine vander_loc(&
!!$         ndeg, ngauss, ninter&
!!$         ix_inter, iy_inter, ix_gauss, iy_gauss,&
!!$         basis1d_onguass_xinter,&
!!$         basis1d_onguass_yinter,&
!!$         vander_loc)
!!$      use Globals
!!$      implicit none
!!$      integer :: ndeg, ngauss, ninter
!!$      integer :: ix_inter, iy_inter
!!$      integer :: ix_gauss, iy_gauss
!!$      real(kind=double) :: basis1d_guassvalues
!!$      real(kind=double) :: vander_loc( (ndeg+1)**2 ) 
!!$      do i = 1, ndeg+1
!!$         do j = 1, ndeg+1
!!$            m = m+1
!!$            vander_loc(mm) = &
!!$                 basis1d_ongauss_xinter(i,ix_gauss) * & 
!!$                 basis1d_ongauss_yinter(j,iy_gauss)
!!$         end do
!!$      end do
!!$    end subroutine vander_loc
!!$
!!$
!!$end function integral_nrm_grad2
!!$
!!$
!!$
!!$subroutine build_grad_quads(nquad, ndeg, &
!!$     coeff, nrm_grad_quad,  basisvalues1d)
!!$  
!!$  call dx2d_coeff( (ndeg+1) **2, coeff, dx_coeff)
!!$  call dy2d_coeff( (ndeg+1) **2, coeff, dy_coeff)
!!$
!!$  do iquad = 1, nquad
!!$     ix_inter = ixinter(iquad, ninter)
!!$     iy_inter = iyinter(iquad, ninter)
!!$     call spec(ix_inter,
!!$
!!$     nrm_grad_quad(iquad) = &
!!$          integral_nrm_grad2(&
!!$          ix_inter,iy_inter,&
!!$          dx_coeff, dy_coeff,&
!!$          basis1d_onguass(:,:,ix_inter),&
!!$          basis1d_onguass(:,:,iy_inter),&
!!$          vander_loc,&
!!$          spec) &
!!$          / area
!!$  end do
!!$
!!$end subroutine build_grad_quads
      
     
  
           
     

