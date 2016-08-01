subroutine structured_grid_analysis
	use common_params
	implicit none
   include 'cgnslib_f.h'
	
   !�J�E���^
	integer:: n				!�P�[�X�����[�v�J�E���^
	integer:: m, id		!���f�ʃ��[�v�J�E���^
	integer:: jp			!���f�ʁE���_���[�v�J�E���^
	integer:: i, j			!�i�q�_���[�v�J�E���^
	integer:: l				!�o�͐����[�v�J�E���^
	integer:: ii			!�o�̓t�@�C�����쐬�p
	!integer:: i, j, k, l, m, n
	
	!�ŒZ�i�q�_�̒T��
	double precision :: ds, dsmin, ss
	
	!�L���|�X�g�l
	double precision :: skp
   
	!�G���[����l
	integer :: ier
	
	!�e�P�[�X�̌v�Z���ʏo�͉񐔁Almax0��BaseCase
	integer :: lmax0, lmax
	
	!��r�i�q�T�C�Y�E�i�q
   integer:: isize0, jsize0, isize, jsize
   double precision, dimension(:,:), allocatable:: xx, yy
	
   !����
   double precision :: tt
	
	!�͏��ϓ��v�Z�t���O
	integer :: ibedcal
    
   !�i�q�����ϐ�
	double precision, dimension(:,:,:), allocatable::  zz0
   double precision, dimension(:,:,:), allocatable::  dep
   double precision, dimension(:,:,:), allocatable::  srf
   double precision, dimension(:,:,:), allocatable::  elv
   double precision, dimension(:,:,:), allocatable::  vx, vy
	double precision, dimension(:,:,:), allocatable:: 	vor
	double precision, dimension(:,:,:), allocatable:: 	fnum
	double precision, dimension(:,:,:), allocatable:: 	sst
	double precision, dimension(:,:,:), allocatable:: 	ech
	double precision, dimension(:,:,:), allocatable:: 	fbe
	double precision, dimension(:,:,:), allocatable:: 	csemin
	double precision, dimension(:,:,:), allocatable:: 	cseave
	double precision, dimension(:,:,:), allocatable:: 	cssave
	double precision, dimension(:,:,:), allocatable:: 	snum
	double precision, dimension(:,:,:), allocatable:: 	blx, bly
	double precision, dimension(:,:,:), allocatable:: 	slc

    
   !���Z�l
   double precision, dimension(:,:,:), allocatable::  delt_dep
   double precision, dimension(:,:,:), allocatable::  delt_srf
   double precision, dimension(:,:,:), allocatable::  delt_elv
   double precision, dimension(:,:,:), allocatable::  delt_vx, delt_vy
	double precision, dimension(:,:,:), allocatable:: 	delt_vor
	double precision, dimension(:,:,:), allocatable:: 	delt_fnum
	double precision, dimension(:,:,:), allocatable:: 	delt_sst
	double precision, dimension(:,:,:), allocatable:: 	delt_ech
	double precision, dimension(:,:,:), allocatable:: 	delt_fbe
	double precision, dimension(:,:,:), allocatable:: 	delt_csemin
	double precision, dimension(:,:,:), allocatable:: 	delt_cseave
	double precision, dimension(:,:,:), allocatable:: 	delt_cssave
	double precision, dimension(:,:,:), allocatable:: 	delt_snum
	double precision, dimension(:,:,:), allocatable:: 	delt_blx, delt_bly
	double precision, dimension(:,:,:), allocatable:: 	delt_slc
	
	!�o�̓t�@�C����
	character(len=22) :: judan_csv
	character(len=33) :: oudan_csv
	
	!work�p
	double precision, dimension(:), allocatable:: wtmp
	integer :: iip, jjp, iip0, jjp0
!================================================================================
	
	!--------------------------------------------------
	!�o�͗p�t�H���_�A�t�@�C�����̐ݒ�
	if(irsd == 1)then
		!�t�H���_�L���m�F
		call system("mkdir "//trim(juoudan_folder)//"\judan")
		call system("del /Q "//trim(juoudan_folder)//"\judan\*.csv")
		call system("mkdir "//trim(juoudan_folder)//"\oudan")
		call system("del /Q "//trim(juoudan_folder)//"\oudan\*.csv")
	
		!�o�̓t�@�C�����f�t�H���g
		judan_csv = 't=00000000.00(sec).csv'
		oudan_csv = 'kp=+0000.0_t=00000000.00(sec).csv'
	end if
	
	!----------------------------------------------------------------------
   !BaseCase�̌v�Z���ʂ���
	! �v�Z���ʏo�͉񐔁Flmax�̐ݒ�  ���@�����o�͉񐔂łȂ��Ɣ�r���Ȃ��@�������A���������ۂ��͔��肵�Ă��Ȃ�
   ! �v�Z���ʏo�͉񐔁Flmax�̐ݒ�  ���@���Ȃ��ق��ɍ��킹��B�������A���������ۂ��͔��肵�Ă��Ȃ�
	
   call cg_iric_read_sol_count_mul_f(idf(0), lmax0, ier)
   
	!�o�͉񐔂̈قȂ��r�P�[�X���Ȃ����G���[�`�F�b�N
   lmax = lmax0
	do n = 1, nmax-1
      call cg_iric_read_sol_count_mul_f(idf(n), lmax0, ier)
      lmax = min(lmax, lmax0)
  !    if(lmax0 /= lmax) then
		!	write(*,*) "***** Warning : Number of Output is mismatch *****"
		!	write(*,*) trim(cname(0))//" : lmax = ", lmax0
		!	write(*,*) trim(cname(n))//" : lmax = ", lmax
		!	goto 997
		!end if
   end do
   !lmax = lmax0
	
	!----------------------------------------------------------------------
   !BaseCase�̌v�Z���ʂ���
	! ��r�̊�b�ƂȂ�i�q:isize, jsize, xx, yy��ݒ肷��@���@cgns�t�@�C���ɏo��
   call cg_iric_gotogridcoord2d_mul_f(idf(0), isize0, jsize0, ier)
	if(ier /= 0) then
		write(*,*) "***** Error : Reading File *****"
		write(*,*) "File = ", trim(fname(0))
		goto 998
	end if
	
	!�i�q���̈قȂ��r�P�[�X���Ȃ����G���[�`�F�b�N
   do n = 1, nmax-1
      call cg_iric_gotogridcoord2d_mul_f(idf(n), isize, jsize, ier)
		if(ier /= 0) then
			write(*,*) "***** Error : Reading File *****"
			write(*,*) "File = ", trim(fname(n))
			goto 998
		end if
		if(isize /= isize0 .or. jsize /= jsize0) then
			write(*,*) "***** Error : Grid Number is Mismatch *****"
			write(*,*) trim(cname(0))//" : isize = ", isize0, ", jsize = ", jsize0
			write(*,*) trim(cname(n))//" : isize = ", isize, ", jsize = ", jsize
			goto 999
		end if
   end do

   !��r����v�Z�Ώۊi�q���m��
   isize = isize0
   jsize = jsize0
	write(*,*) "Grid isize", isize
	write(*,*) "Grid jsize", jsize
	
   allocate(xx(1:isize, 1:jsize), yy(1:isize, 1:jsize))
   call cg_iric_getgridcoord2d_mul_f(idf(0), xx, yy, ier)
   call cg_iric_writegridcoord2d_mul_f(idf0, isize, jsize, xx, yy, ier)
	
	
	!----------------------------------------------------------------------
	!�c�f�}�p�@
	!���@���S�_�A����сA�e���f���_�ɍł��߂��i�q�_id��ݒ�
	if(irsd == 1)then
	
		!�͓����S�_�ɍł��߂��i�q�_�ԍ�
		do m = 1, mmax
			dsmin = 9999
			judan(m)%ip = -1
			judan(m)%jp = -1
			
			do j = 1, jsize
				do i = 1, isize
					ds = sqrt((xx(i,j) - judan(m)%xx)**2. + (yy(i,j) - judan(m)%yy)**2.)
					if(ds < dsmin)then
						dsmin = ds
						judan(m)%ip = i
						judan(m)%jp = j
					end if
				end do
			end do
		end do
	end if
	
	!----------------------------------------------------------------------
	!���f�}�p�@
	!���@���S�_�A����сA�e���f���_�ɍł��߂��i�q�_id��ݒ�
	if(irsd == 1)then
	
		!�e���f���_�ɍł��߂��i�q�_�ԍ�
		do m = 1, mmax  !�f�ʃ��[�v
			do jp = 1, odan(m)%jpmax  !���_���[�v
				dsmin = 9999
				odan(m)%ip(jp) = -1
				odan(m)%jp(jp) = -1
				do j = 1, jsize
					do i = 1, isize
						ds = sqrt((xx(i,j) - odan(m)%xx(jp))**2. + (yy(i,j) - odan(m)%yy(jp))**2.)
						if(ds < dsmin)then
							dsmin = ds
							odan(m)%ip(jp) = i
							odan(m)%jp(jp) = j
						end if
					end do
				end do
			end do
		end do
	end if
   
	!�����l�̏o��
	call rescomp(lmax, isize, jsize, xx, yy)
		
	return
	
999	write(*,*) "CAT2D can be applied only to the same grid size condition."
		stop
		
998	write(*,*) "Structured is specified as Grid type condition."
		write(*,*) "But above CGNS file does not have structured Grid."
		write(*,*) "Please check your calculation condition."
		stop

997	write(*,*) "CAT2D can be applied only to the same number of output condition."
		stop
		
		
end subroutine