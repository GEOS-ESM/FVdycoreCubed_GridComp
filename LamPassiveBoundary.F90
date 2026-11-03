! Isolated manufactured-passive contract probe for the dry standalone host.
module LamPassiveBoundaryMod
  use ESMF
  use netcdf
  use FV_StateMod, only: FV_Atm
  use fv_lam_geometry_mod, only: lam_nc_check
  use fv_tracer2d_mod, only: set_buffer_probe
  implicit none
  private
  public :: lam_passive_initial, lam_passive_before_adv, lam_passive_checkpoint
  integer, save :: cached_lower=-1, cached_upper=-1
  real(kind=8), allocatable, save :: low(:,:,:,:), high(:,:,:,:)
contains
  subroutine read_field(file, name, values)
    character(*), intent(in) :: file, name
    real(kind=8), intent(out) :: values(:,:,:)
    integer :: nc, vid, st, var_rank, dims(3), n, length
    character(32) :: units
    call lam_nc_check(nf90_open(file,nf90_nowrite,nc))
    call lam_nc_check(nf90_inq_varid(nc,name,vid))
    call lam_nc_check(nf90_inquire_variable(nc,vid,ndims=var_rank,dimids=dims))
    if(var_rank/=3) error stop 'LAM passive file rank'
    call lam_nc_check(nf90_inquire_dimension(nc,dims(1),len=length))
    if(length/=FV_Atm(1)%npx-1+8) error stop 'LAM passive x extent'
    call lam_nc_check(nf90_inquire_dimension(nc,dims(2),len=length))
    if(length/=FV_Atm(1)%npy-1+8) error stop 'LAM passive y extent'
    call lam_nc_check(nf90_inquire_dimension(nc,dims(3),len=length))
    if(length/=FV_Atm(1)%npz) error stop 'LAM passive levels'
    call lam_nc_check(nf90_get_att(nc,vid,'units',units))
    if(trim(units)/='kg kg-1') error stop 'LAM passive units'
    call lam_nc_check(nf90_get_var(nc,vid,values, &
         start=[FV_Atm(1)%bd%is+4,FV_Atm(1)%bd%js+4,1]))
    if(any(values<0.).or.any(values/=values)) error stop 'LAM passive input invalid'
    call lam_nc_check(nf90_close(nc))
  end subroutine

  subroutine lam_passive_checkpoint(bundle)
    type(ESMF_FieldBundle), intent(in) :: bundle
    type(ESMF_VM) :: vm
    type(ESMF_Field) :: field
    real(kind=4), pointer :: q(:,:,:)
    character(ESMF_MAXSTR) :: name
    character(256) :: file
    integer :: status,n,nq,rank,unit,seconds,shape3(3)
    call ESMF_FieldBundleGet(bundle,fieldCount=nq,rc=status)
    if(status/=ESMF_SUCCESS) error stop 'LAM passive checkpoint bundle'
    if(nq<=1) return
    seconds=nint(FV_Atm(1)%neststruct%lam_elapsed)
    if(mod(seconds,1800)/=0) return
    call ESMF_VMGetGlobal(vm,rc=status)
    if(status/=ESMF_SUCCESS) error stop 'LAM passive checkpoint VM'
    call ESMF_VMGet(vm,localPet=rank,rc=status)
    if(status/=ESMF_SUCCESS) error stop 'LAM passive checkpoint rank'
    write(file,'("lam-passive-",I6.6,"-p",I4.4,".bin")')seconds,rank
    open(newunit=unit,file=trim(file),status='new',access='stream',form='unformatted',action='write',iostat=status)
    if(status/=0) error stop 'LAM passive checkpoint file exists or unwritable'
    write(unit) 440314,seconds,nq
    do n=1,nq
       call ESMF_FieldBundleGet(bundle,fieldIndex=n,field=field,rc=status)
       if(status/=ESMF_SUCCESS) error stop 'LAM passive checkpoint field'
       call ESMF_FieldGet(field,name=name,rc=status)
       if(status/=ESMF_SUCCESS) error stop 'LAM passive checkpoint name'
       call ESMF_FieldGet(field,farrayPtr=q,rc=status)
       if(status/=ESMF_SUCCESS) error stop 'LAM passive checkpoint pointer'
       shape3=shape(q)
       write(unit) name(1:4),shape3,q
    enddo
    close(unit)
    print *, 'LAM_PASSIVE_CHECKPOINT',seconds,rank,nq,trim(file)
  end subroutine

  subroutine lam_passive_initial(bundle)
    type(ESMF_FieldBundle), intent(in) :: bundle
    type(ESMF_Field) :: field
    real(kind=4), pointer :: q(:,:,:)
    real(kind=8), allocatable :: tmp(:,:,:)
    character(ESMF_MAXSTR) :: name
    character(256) :: file
    integer :: status,n,nq
    call ESMF_FieldBundleGet(bundle,fieldCount=nq,rc=status)
    if(status/=ESMF_SUCCESS) error stop 'LAM passive initial bundle'
    if(nq<=1) return
    file=trim(FV_Atm(1)%neststruct%lam_directory)//'/lam.initial.nc'
    do n=2,nq
       call ESMF_FieldBundleGet(bundle,fieldIndex=n,field=field,rc=status)
       if(status/=ESMF_SUCCESS) error stop 'LAM passive initial field'
       call ESMF_FieldGet(field,name=name,rc=status)
       if(status/=ESMF_SUCCESS) error stop 'LAM passive initial name'
       call ESMF_FieldGet(field,farrayPtr=q,rc=status)
       if(status/=ESMF_SUCCESS) error stop 'LAM passive initial pointer'
       allocate(tmp(size(q,1),size(q,2),size(q,3)))
       write(name,'("q",I3.3)') n-1
       call read_field(trim(file),trim(name),tmp)
       q=real(tmp,kind=4)
       print *, 'LAM_PASSIVE_INITIAL',n,trim(name),minval(q),maxval(q)
       deallocate(tmp)
    enddo
  end subroutine

  subroutine lam_passive_before_adv(bundle,ple0,clock)
    type(ESMF_FieldBundle), intent(in) :: bundle
    real(kind=8), intent(in) :: ple0(:,:,:)
    type(ESMF_Clock), intent(in) :: clock
    type(ESMF_Time) :: now,origin
    type(ESMF_TimeInterval) :: gap
    type(ESMF_Field) :: field
    real(kind=4), pointer :: q(:,:,:)
    real(kind=8) :: elapsed, alpha, imposed, value, area, dp
    integer :: status,n,nq,i,j,k,gi,gj,lower,upper,interval
    character(4) :: var
    character(256) :: file
    call ESMF_FieldBundleGet(bundle,fieldCount=nq,rc=status)
    if(status/=ESMF_SUCCESS) error stop 'LAM passive ADV bundle'
    if(nq<=1) return
    associate(o=>FV_Atm(1)%neststruct%lam_origin)
      call ESMF_TimeSet(origin,yy=o(1),mm=o(2),dd=o(3),h=o(4),m=o(5),s=o(6),rc=status)
    end associate
    if(status/=ESMF_SUCCESS) error stop 'LAM passive clock origin'
    call ESMF_ClockGet(clock,currTime=now,rc=status)
    if(status/=ESMF_SUCCESS) error stop 'LAM passive clock now'
    gap=now-origin
    call ESMF_TimeIntervalGet(gap,s_r8=elapsed,rc=status)
    if(status/=ESMF_SUCCESS) error stop 'LAM passive clock seconds'
    print *, 'GBUF_CLOCK',elapsed,FV_Atm(1)%neststruct%lam_elapsed
    interval=FV_Atm(1)%neststruct%lam_interval
    lower=floor(elapsed/(3600.d0*interval))*interval
    upper=lower+interval
    if(elapsed<lower*3600.d0.or.elapsed>upper*3600.d0) error stop 'LAM passive stage time'
    alpha=(elapsed-lower*3600.d0)/(3600.d0*interval)
    if(lower/=cached_lower.or.upper/=cached_upper) then
       call ESMF_FieldBundleGet(bundle,fieldIndex=2,field=field,rc=status)
       if(status/=ESMF_SUCCESS) error stop 'LAM passive first field'
       call ESMF_FieldGet(field,farrayPtr=q,rc=status)
       if(status/=ESMF_SUCCESS) error stop 'LAM passive first pointer'
       if(allocated(low)) deallocate(low,high)
       allocate(low(size(q,1),size(q,2),size(q,3),nq-1))
       allocate(high(size(q,1),size(q,2),size(q,3),nq-1))
       do n=2,nq
          write(var,'("q",I3.3)')n-1
          write(file,'(A,"/lam.",I3.3,".nc")')trim(FV_Atm(1)%neststruct%lam_directory),lower
          call read_field(trim(file),var,low(:,:,:,n-1))
          write(file,'(A,"/lam.",I3.3,".nc")')trim(FV_Atm(1)%neststruct%lam_directory),upper
          call read_field(trim(file),var,high(:,:,:,n-1))
       enddo
       cached_lower=lower;cached_upper=upper
       print *, 'LAM_PASSIVE_BRACKET',elapsed,lower*3600,upper*3600
    endif
    call set_buffer_probe(low,high,real(lower*3600,kind=8), &
         real(upper*3600,kind=8),elapsed)
    if(any(shape(ple0) /= [size(low,1),size(low,2),size(low,3)+1])) error stop 'LAM passive PLE0 extent'
    do n=2,nq
       call ESMF_FieldBundleGet(bundle,fieldIndex=n,field=field,rc=status)
       if(status/=ESMF_SUCCESS) error stop 'LAM passive boundary field'
       call ESMF_FieldGet(field,farrayPtr=q,rc=status)
       if(status/=ESMF_SUCCESS) error stop 'LAM passive boundary pointer'
       imposed=0.d0
       do k=1,size(q,3)
       do j=1,size(q,2)
       do i=1,size(q,1)
          gi=FV_Atm(1)%bd%is+i-1;gj=FV_Atm(1)%bd%js+j-1
          if(gi<=3.or.gi>=FV_Atm(1)%npx-3.or.gj<=3.or.gj>=FV_Atm(1)%npy-3)then
             value=(1.d0-alpha)*low(i,j,k,n-1)+alpha*high(i,j,k,n-1)
             area=FV_Atm(1)%gridstruct%area_64(gi,gj)
             dp=ple0(i,j,k+1)-ple0(i,j,k)
             imposed=imposed+(real(value,kind=4)-q(i,j,k))*dp*area/9.80665d0
             q(i,j,k)=real(value,kind=4)
          endif
       enddo
       enddo
       enddo
       print *, 'LAM_PASSIVE_BOUNDARY',n,elapsed,alpha,imposed,minval(q),maxval(q)
    enddo
  end subroutine
end module
