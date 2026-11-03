!============ GEOS-LAM Update Start =============
! Construct the one-tile bounded MAPL grid and staggering.
module LamGridMod
 use ESMF
 use MAPL
 use MAPL_XYGridFactoryMod,only:XYGridFactory
 use MAPL_KeywordEnforcerMod,only:KeywordEnforcer
 use FV_StateMod,only:FV_Atm
 implicit none
 type,extends(XYGridFactory)::LamGridFactory
    type(ESMF_Grid)::native_grid
 contains
    procedure::make_new_grid=>lam_factory_make
 end type
contains
 function lam_factory_make(this,unusable,rc) result(grid)
 class(LamGridFactory),intent(in)::this
 class(KeywordEnforcer),optional,intent(in)::unusable
 integer,optional,intent(out)::rc
 type(ESMF_Grid)::grid
 grid=this%native_grid
 if(present(rc))rc=ESMF_SUCCESS
 end function

 subroutine create_lam_grid(gc,rc)
 type(ESMF_GridComp),intent(inout)::gc
 integer,intent(out)::rc
 type(ESMF_Grid)::grid
 type(LamGridFactory)::factory
 type(MAPL_MetaComp),pointer::maplobj
 real(ESMF_KIND_R8),pointer::x(:,:),y(:,:)
 integer::nx,ny,i,j,si,sj,status
 integer,allocatable::ims(:),jms(:)
 real(ESMF_KIND_R8),allocatable::cornx(:),corny(:)
 nx=FV_Atm(1)%layout(1);ny=FV_Atm(1)%layout(2)
 if(mod(FV_Atm(1)%npx-1,nx)/=0.or.mod(FV_Atm(1)%npy-1,ny)/=0) error stop 'LAM requires divisible MAPL/FMS layout'
 allocate(ims(nx),jms(ny));ims=(FV_Atm(1)%npx-1)/nx;jms=(FV_Atm(1)%npy-1)/ny
 grid=ESMF_GridCreateNoPeriDim(name='NASA-Standalone-LAM',countsPerDEDim1=ims,countsPerDEDim2=jms, &
 indexFlag=ESMF_INDEX_DELOCAL,gridEdgeLWidth=[0,0],gridEdgeUWidth=[1,1], &
 coordDep1=[1,2],coordDep2=[1,2],coordSys=ESMF_COORDSYS_SPH_RAD,coordTypeKind=ESMF_TYPEKIND_R8,rc=status)
 if(status/=ESMF_SUCCESS) error stop 'LAM grid create'
 call ESMF_GridAddCoord(grid,staggerEdgeLWidth=[0,0],staggerEdgeUWidth=[0,0],rc=status)
 if(status/=ESMF_SUCCESS)error stop 'LAM previous ESMF operation failed'
 call ESMF_GridGetCoord(grid,coordDim=1,farrayPtr=x,rc=status)
 if(status/=ESMF_SUCCESS)error stop 'LAM previous ESMF operation failed'
 call ESMF_GridGetCoord(grid,coordDim=2,farrayPtr=y,rc=status)
 si=FV_Atm(1)%bd%is;sj=FV_Atm(1)%bd%js
 if(size(x,1)/=FV_Atm(1)%bd%ie-si+1.or.size(x,2)/=FV_Atm(1)%bd%je-sj+1) error stop 'LAM ESMF/FMS domain mismatch'
 if(status/=ESMF_SUCCESS)error stop 'LAM coordinates unavailable'
 do j=1,size(x,2)
 do i=1,size(x,1)
 x(i,j)=FV_Atm(1)%gridstruct%agrid(si+i-1,sj+j-1,1)
 y(i,j)=FV_Atm(1)%gridstruct%agrid(si+i-1,sj+j-1,2)
 enddo
 enddo
 allocate(cornx((size(x,1)+1)*(size(x,2)+1)),corny((size(x,1)+1)*(size(x,2)+1)))
 do j=1,size(x,2)+1
 do i=1,size(x,1)+1
 cornx((j-1)*(size(x,1)+1)+i)=FV_Atm(1)%gridstruct%grid(si+i-1,sj+j-1,1)
 corny((j-1)*(size(x,1)+1)+i)=FV_Atm(1)%gridstruct%grid(si+i-1,sj+j-1,2)
 enddo
 enddo
 call ESMF_AttributeSet(grid,'GridType','XY',rc=status)
 call ESMF_AttributeSet(grid,'GRID_LM',FV_Atm(1)%npz,rc=status)
 call ESMF_AttributeSet(grid,'GridCornerLons:',itemCount=size(cornx),valueList=cornx,rc=status)
 call ESMF_AttributeSet(grid,'GridCornerLats:',itemCount=size(corny),valueList=corny,rc=status)
 call ESMF_GridAddCoord(grid,staggerloc=ESMF_STAGGERLOC_CORNER,staggerEdgeLWidth=[0,0],staggerEdgeUWidth=[1,1],rc=status)
 if(status/=ESMF_SUCCESS)error stop 'LAM previous ESMF operation failed'
 call ESMF_GridGetCoord(grid,coordDim=1,staggerloc=ESMF_STAGGERLOC_CORNER,farrayPtr=x,rc=status)
 if(status/=ESMF_SUCCESS)error stop 'LAM previous ESMF operation failed'
 call ESMF_GridGetCoord(grid,coordDim=2,staggerloc=ESMF_STAGGERLOC_CORNER,farrayPtr=y,rc=status)
 if(status/=ESMF_SUCCESS)error stop 'LAM corners unavailable'
 if(size(x,1)<ims(1).or.size(x,2)<jms(1))error stop 'LAM corner extent'
 do j=1,size(x,2)
 do i=1,size(x,1)
 x(i,j)=FV_Atm(1)%gridstruct%grid(si+i-1,sj+j-1,1)
 y(i,j)=FV_Atm(1)%gridstruct%grid(si+i-1,sj+j-1,2)
 enddo
 enddo
 factory%XYGridFactory=XYGridFactory(grid_file_name='INPUT/lam-mapl-grid.nc',grid_name='NASA-Standalone-LAM', &
 im_world=FV_Atm(1)%npx-1,jm_world=FV_Atm(1)%npy-1,lm=FV_Atm(1)%npz,nx=nx,ny=ny,rc=status)
 if(status/=ESMF_SUCCESS) error stop 'LAM MAPL XY metadata init'
 factory%native_grid=grid
 grid=grid_manager%make_grid(factory,rc=status)
 if(status/=ESMF_SUCCESS) error stop 'LAM MAPL factory registration'
 call MAPL_GetObjectFromGC(gc,maplobj,rc=status)
 call maplobj%grid%set(grid,rc=status)
 call ESMF_GridCompSet(gc,grid=grid,rc=rc)
 end subroutine
end module
!============ GEOS-LAM Update End =============
