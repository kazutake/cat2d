subroutine unstructured_grid_analysis
	use common_params
	use grid
	
	implicit none
   include 'cgnslib_f.h'
	
   !�J�E���^
	integer:: n				!�P�[�X�����[�v�J�E���^
	integer:: m, id		!���f�ʃ��[�v�J�E���^
	integer:: jp			!���f�ʁE���_���[�v�J�E���^
	integer:: np, np0			!�i�q�_���[�v�J�E���^
	integer:: l				!�o�͐����[�v�J�E���^
	integer:: ii			!�o�̓t�@�C�����쐬�p
	
	!�ŒZ�i�q�_�̒T��
	double precision :: ds, dsmin, ss
	
	!�L���|�X�g�l
	double precision :: skp

	!�G���[����l
	integer :: ier
	
	!�e�P�[�X�̌v�Z���ʏo�͉񐔁Almax0��BaseCase
	integer :: lmax0, lmax
	
	!�i�q�R�s�[�p�ϐ�
	integer :: baseId, indexId, zoneId, gridId
	character(len=strmax) :: zonename
	integer :: zoneCopied, gridCopied
	integer :: n_pts, n_elems, nsections
	
	!�i�q�T�C�Y�E�i�q�E�v�f
   integer, dimension(3,3) :: isize, isize0
   double precision, dimension(:), allocatable:: xx, yy
	integer, dimension(:), allocatable :: element, parentdata
      
   !����
   double precision :: tt
	
	!�͏��ϓ��v�Z�t���O
	integer :: ibedcal
	
   !�i�q�����ϐ�
   double precision, dimension(:,:), allocatable:: dep
	double precision, dimension(:,:), allocatable:: zz0
   double precision, dimension(:,:), allocatable:: srf
   double precision, dimension(:,:), allocatable:: frv
   double precision, dimension(:,:), allocatable:: flx
   double precision, dimension(:,:), allocatable:: vx, vy
	double precision, dimension(:,:), allocatable:: elv
	double precision, dimension(:,:), allocatable:: ech
	double precision, dimension(:,:), allocatable:: blx, bly
	double precision, dimension(:,:), allocatable:: cgd
	double precision, dimension(:,:), allocatable:: vbx, vby
	double precision, dimension(:,:), allocatable:: slc
	
   !���Z�l
	double precision, dimension(:,:), allocatable:: delt_dep
	double precision, dimension(:,:), allocatable:: delt_zz0
   double precision, dimension(:,:), allocatable:: delt_srf
   double precision, dimension(:,:), allocatable:: delt_frv
   double precision, dimension(:,:), allocatable:: delt_flx
   double precision, dimension(:,:), allocatable:: delt_vx, delt_vy
	double precision, dimension(:,:), allocatable:: delt_elv
	double precision, dimension(:,:), allocatable:: delt_ech
	double precision, dimension(:,:), allocatable:: delt_blx, delt_bly
	double precision, dimension(:,:), allocatable:: delt_cgd
	double precision, dimension(:,:), allocatable:: delt_vbx, delt_vby
	double precision, dimension(:,:), allocatable:: delt_slc
	
	!�o�̓t�@�C����
	character(len=22) :: judan_csv
	character(len=33) :: oudan_csv
	
	!work�p�ϐ�
	double precision, dimension(:), allocatable:: wtmp
	!integer :: nnp, nnp0
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
	
	!--------------------------------------------------
   !�v�Z�i�q�����݂�cgns�t�@�C���ɃR�s�[
	call cg_nbases_f(idf(0), baseId, ier)  !Get number of CGNS base nodes in file
	call cg_nzones_f(idf(0), baseId, zoneId, ier)  !Get number of zones in base
	!baseId = 1
	!zoneId = 1
	gridId = 1
	
	!�R�s�[
	call copyZoneIfNeeded(idf(0), idf0, baseId, zoneId, zoneCopied)
	call copyGrid(idf(0), idf0, baseId, zoneId, gridId, gridId, gridCopied)
	call copyAllElements(idf(0), idf0, baseId, zoneId)
	call cg_iric_initgrid_mul_f(idf0, zoneId, ier)	!�R�s�[�i�q��iriclib �ŗ��p�ł���悤������
	
	!--------------------------------------------------
	!����̊i�q�����p����Ă��邱�Ƃ��m�F
	call cg_zone_read_f(idf0, baseId, zoneId, zonename, isize0, ier)
	if(ier /= 0) then
		write(*,*) "***** Error : Reading File *****"
		write(*,*) "File = ", trim(fname(n))
		goto 998
	end if
	do n = 0, nmax-1
		call cg_zone_read_f(idf(n), baseId, zoneId, zonename, isize, ier)
		if(ier /= 0) then
			write(*,*) "***** Error : Reading File *****"
			write(*,*) "File = ", trim(fname(n))
			goto 998
		end if
		if(isize(1,1) /= isize0(1,1) .or. isize(2,1) /= isize0(2,1)) then
			write(*,*) "***** Error : Grid Number is Mismatch *****"
			write(*,*) trim(cname(0))//" : Pts = ", isize0(1,1), ", Elems = ", isize0(2,1)
			write(*,*) trim(cname(n))//" : Pts = ", isize(1,1),  ", Elems = ", isize(2,1)
			goto 999
		end if
	end do
	
	!--------------------------------------------------
	!�i�q����
	n_pts = isize0(1,1)
	n_elems = isize0(2,1)
	write(*,*) "Number of Ponts    = ", n_pts
	write(*,*) "Number of Elements = ", n_elems
	
	!�i�q���W��ǂݍ���
	allocate(xx(1:n_pts), yy(1:n_pts))
	call cg_coord_read_f(idf0, baseId, zoneId, 'CoordinateX', RealDouble, 1, n_pts, xx, ier)
	call cg_coord_read_f(idf0, baseId, zoneId, 'CoordinateY', RealDouble, 1, n_pts, yy, ier)
	
	!----------------------------------------------------------------------
   !�v�Z���ʏo�͉񐔂̊m�F  
	! ���@�����o�͉񐔁A�������łȂ��Ɣ�r���Ȃ��@�������A���������ۂ��͔��肵�Ȃ�
   call cg_iric_read_sol_count_mul_f(idf(0), lmax0, ier)
   do n = 1, nmax-1
      call cg_iric_read_sol_count_mul_f(idf(n), lmax, ier)
      if(lmax0 /= lmax) then
			write(*,*) "***** Error : Number of Output is mismatch *****"
			write(*,*) trim(cname(0))//" : lmax = ", lmax0
			write(*,*) trim(cname(n))//" : lmax = ", lmax
			goto 997
		end if
   end do
   lmax = lmax0
	
	!--------------------------------------------------
	!�c�f�}�p�@���@���S�_�A����сA�e���f���_�ɍł��߂��i�q�_id��ݒ�
	if(irsd == 1)then
	
		!�͓����S�_�ɍł��߂��i�q�_�ԍ�
		do m = 1, mmax
			dsmin = 9999
			judan(m)%np = -1
			
			do np = 1, n_pts
				ds = sqrt((xx(np) - judan(m)%xx)**2. + (yy(np) - judan(m)%yy)**2.)
				if(ds < dsmin)then
					dsmin = ds
					judan(m)%np = np
				end if
		
			end do
		end do
	
		!�e���f���_�ɍł��߂��i�q�_�ԍ�
		do m = 1, mmax
			do jp = 1, odan(m)%jpmax
				dsmin = 9999
				odan(m)%np(jp) = -1
				do np = 1, n_pts
					ds = sqrt((xx(np) - odan(m)%xx(jp))**2. + (yy(np) - odan(m)%yy(jp))**2.)
					if(ds < dsmin)then
						dsmin = ds
						odan(m)%np(jp) = np
					end if
				end do
			end do
		end do
	end if
	
   !--------------------------------------------------
   !�ϐ��A���P�[�V����
   allocate(dep(0:nmax-1,n_pts))
   allocate(srf(0:nmax-1,n_pts))
   allocate(zz0(0:nmax-1,n_pts))
   allocate(flx(0:nmax-1,n_pts))
	allocate(frv(0:nmax-1,n_pts))
   allocate(vx(0:nmax-1,n_pts), vy(0:nmax,n_pts))
	
	allocate(elv(0:nmax-1,n_pts))
	allocate(ech(0:nmax-1,n_pts))
	allocate(blx(0:nmax-1,n_pts), bly(0:nmax-1,n_pts))
	allocate(cgd(0:nmax-1,n_pts))
	allocate(vbx(0:nmax-1,n_pts), vby(0:nmax-1,n_pts))
	allocate(slc(0:nmax-1,n_pts))
	
	allocate(delt_dep(0:nmax-1,n_pts))
	allocate(delt_zz0(0:nmax-1,n_pts))
   allocate(delt_srf(0:nmax-1,n_pts))
   allocate(delt_frv(0:nmax-1,n_pts))
   allocate(delt_flx(0:nmax-1,n_pts))
   allocate(delt_vx(0:nmax-1,n_pts), delt_vy(0:nmax-1,n_pts))
	allocate(delt_elv(0:nmax-1,n_pts))
	allocate(delt_ech(0:nmax-1,n_pts))
	allocate(delt_blx(0:nmax-1,n_pts), delt_bly(0:nmax-1,n_pts))
	allocate(delt_cgd(0:nmax-1,n_pts))
	allocate(delt_vbx(0:nmax-1,n_pts), delt_vby(0:nmax-1,n_pts))
	allocate(delt_slc(0:nmax-1,n_pts))
	
	allocate(wtmp(n_pts))
		
   !----------------------------------------------------------------------
	!�v�Z���ʓǍ��J�n
   do l = 1, lmax
	 
		!--------------------------------------------------
      !�f�[�^�Ǎ��@�����ʂ�Ǎ���
		!--------------------------------------------------
		call cg_iric_read_sol_time_mul_f(idf(0), l, tt, ier)		!��{�f�[�^��base case����ǂݍ���
      do n = 0, nmax-1	
         call cg_iric_read_sol_real_mul_f(idf(n), l, "Initial Elevation (m)", zz0(n,:), ier)
			call cg_iric_read_sol_real_mul_f(idf(n), l, "Water surface elevation (m)", srf(n,:), ier)
			call cg_iric_read_sol_real_mul_f(idf(n), l, "Depth (m)", dep(n,:), ier)
			call cg_iric_read_sol_real_mul_f(idf(n), l, "Flux (m**2_s)", flx(n,:), ier)
			call cg_iric_read_sol_real_mul_f(idf(n), l, "Friction velocity (m_s)", frv(n,:), ier)
         call cg_iric_read_sol_real_mul_f(idf(n), l, "Velocity (m_s)_X", vx(n,:), ier)
         call cg_iric_read_sol_real_mul_f(idf(n), l, "Velocity (m_s)_Y", vy(n,:), ier)
			
			!�͏��ϓ��v�Z����ꍇ
			call cg_iric_read_integer_mul_f(idf(n), "ICON_13", ibedcal, ier)  !Mflow�̂�
			if(ibedcal == 1)then
				call cg_iric_read_sol_real_mul_f(idf(n), l, "Elevation (m)",               elv(n,:), ier)
				call cg_iric_read_sol_real_mul_f(idf(n), l, "Elevation change (m)",        ech(n,:), ier)
				call cg_iric_read_sol_real_mul_f(idf(n), l, "Bedload (mm**2_s)_X",         blx(n,:), ier)
				call cg_iric_read_sol_real_mul_f(idf(n), l, "Bedload (mm**2_s)_Y",         bly(n,:), ier)
				call cg_iric_read_sol_real_mul_f(idf(n), l, "Central grain diameter (mm)", cgd(n,:), ier)
				call cg_iric_read_sol_real_mul_f(idf(n), l, "Vel(bottom) (m_s)_X",         vbx(n,:), ier)
				call cg_iric_read_sol_real_mul_f(idf(n), l, "Vel(bottom) (m_s)_Y",         vby(n,:), ier)
				call cg_iric_read_sol_real_mul_f(idf(n), l, "Suspended load",              slc(n,:), ier)
			end if
		end do
		
		!--------------------------------------------------
      !�����l�̎Z��
		!--------------------------------------------------
      do n = 1, nmax-1
			do np = 1, n_pts
				delt_dep(n,np)  = dep(n,np) - dep(0,np) 
				delt_zz0(n,np) = zz0(n,np) - zz0(0,np)
				delt_srf(n,np) = srf(n,np) - srf(0,np)
				delt_frv(n,np) = frv(n,np) - frv(0,np)
				delt_flx(n,np) = flx(n,np) - flx(0,np)
				delt_vx(n,np) = vx(n,np) - vx(0,np)
				delt_vy(n,np) = vy(n,np) - vy(0,np)
				delt_elv(n,np) = elv(n,np) - elv(0,np)
				delt_ech(n,np) = ech(n,np) - ech(0,np)
				delt_blx(n,np) = blx(n,np) - blx(0,np)
				delt_bly(n,np) = bly(n,np) - bly(0,np)
				delt_cgd(n,np) = cgd(n,np) - cgd(0,np)
				delt_vbx(n,np) = vbx(n,np) - vbx(0,np)
				delt_vby(n,np) = vby(n,np) - vby(0,np)
				delt_slc(n,np) = slc(n,np) - slc(0,np)
			end do
	   end do
		
		!--------------------------------------------------
		!�o��
		!--------------------------------------------------
		call cg_iric_write_sol_time_mul_f(idf0, tt, ier)
      write(*,*) "time = ", tt, "[sec]"
		
		!��r�ł���悤�@���̂܂܂��ׂĂ��o��
		do n = 0, nmax-1
			call cg_iric_write_sol_real_mul_f(idf0, trim(cname(n))//"_InitialElevation", zz0(n,:), ier)
			call cg_iric_write_sol_real_mul_f(idf0, trim(cname(n))//"_WaterSurfaceElevation", srf(n,:), ier)
			call cg_iric_write_sol_real_mul_f(idf0, trim(cname(n))//"_Depth", dep(n,:), ier)
			call cg_iric_write_sol_real_mul_f(idf0, trim(cname(n))//"_Flux", flx(n,:), ier)
			call cg_iric_write_sol_real_mul_f(idf0, trim(cname(n))//"_FrictionVelocity", frv(n,:), ier)
         call cg_iric_write_sol_real_mul_f(idf0, trim(cname(n))//"_VelocityX", vx(n,:), ier)
         call cg_iric_write_sol_real_mul_f(idf0, trim(cname(n))//"_VelocityY", vy(n,:), ier)
			
			if(ibedcal == 1)then
				call cg_iric_write_sol_real_mul_f(idf0, trim(cname(n))//"_Elevation",               elv(n,:), ier)
				call cg_iric_write_sol_real_mul_f(idf0, trim(cname(n))//"_ElevationChange",        ech(n,:), ier)
				call cg_iric_write_sol_real_mul_f(idf0, trim(cname(n))//"_Bedload_X",         blx(n,:), ier)
				call cg_iric_write_sol_real_mul_f(idf0, trim(cname(n))//"_Bedload_Y",         bly(n,:), ier)
				call cg_iric_write_sol_real_mul_f(idf0, trim(cname(n))//"_CentralGrainDiameter", cgd(n,:), ier)
				call cg_iric_write_sol_real_mul_f(idf0, trim(cname(n))//"_Vel(bottom)X",         vbx(n,:), ier)
				call cg_iric_write_sol_real_mul_f(idf0, trim(cname(n))//"_Vel(bottom)Y",         vby(n,:), ier)
				call cg_iric_write_sol_real_mul_f(idf0, trim(cname(n))//"_Suspendedload",              slc(n,:), ier)
			end if

      end do
		
		!�����l���o��
		do n = 1, nmax-1
			call cg_iric_write_sol_real_mul_f(idf0, trim(cname(n))//"-"//trim(cname(0))//"_InitialElevation", delt_zz0(n,:), ier)	
			call cg_iric_write_sol_real_mul_f(idf0, trim(cname(n))//"-"//trim(cname(0))//"_WaterSurfaceElevation", delt_srf(n,:), ier)
			call cg_iric_write_sol_real_mul_f(idf0, trim(cname(n))//"-"//trim(cname(0))//"_Depth", delt_dep(n,:), ier)
			call cg_iric_write_sol_real_mul_f(idf0, trim(cname(n))//"-"//trim(cname(0))//"_Flux", delt_flx(n,:), ier)
			call cg_iric_write_sol_real_mul_f(idf0, trim(cname(n))//"-"//trim(cname(0))//"_FrictionVelocity", delt_frv(n,:), ier)
			call cg_iric_write_sol_real_mul_f(idf0, trim(cname(n))//"-"//trim(cname(0))//"_VelocityX", delt_vx(n,:), ier)
			call cg_iric_write_sol_real_mul_f(idf0, trim(cname(n))//"-"//trim(cname(0))//"_VelocityY", delt_vy(n,:), ier)
			
			call cg_iric_write_sol_real_mul_f(idf0, trim(cname(n))//"-"//trim(cname(0))//"_Elevation", delt_elv(n,:), ier)
			call cg_iric_write_sol_real_mul_f(idf0, trim(cname(n))//"-"//trim(cname(0))//"_ElevationChange", delt_ech(n,:), ier)
			call cg_iric_write_sol_real_mul_f(idf0, trim(cname(n))//"-"//trim(cname(0))//"_BedloadX", delt_blx(n,:), ier)
			call cg_iric_write_sol_real_mul_f(idf0, trim(cname(n))//"-"//trim(cname(0))//"_BedloadY", delt_bly(n,:), ier)
			call cg_iric_write_sol_real_mul_f(idf0, trim(cname(n))//"-"//trim(cname(0))//"_CentralGrainDiameter", delt_cgd(n,:), ier)
			call cg_iric_write_sol_real_mul_f(idf0, trim(cname(n))//"-"//trim(cname(0))//"_Vel(bottom)X", delt_vbx(n,:), ier)
			call cg_iric_write_sol_real_mul_f(idf0, trim(cname(n))//"-"//trim(cname(0))//"_Vel(bottom)Y", delt_vby(n,:), ier)
			call cg_iric_write_sol_real_mul_f(idf0, trim(cname(n))//"-"//trim(cname(0))//"_Suspendedload", delt_slc(n,:), ier)
		end do
		
		if(irsd == 1) then
		
			!======================================================================
			!�ȉ��A�c���f�}�쐬�p�f�[�^��csv�`���ŏo��
			!======================================================================
		
			!----------------------------------------------------------------------
			!�c�f�}�p�f�[�^�̏o��
			!judan_csv = 't=00000000.00(sec).csv'
			write(judan_csv(3:13),'(f10.2)') tt
			do ii=3,13
				if(judan_csv(ii:ii) == ' ') judan_csv(ii:ii)='0'
			end do
		
			open(10, file=trim(juoudan_folder)//"\judan\"//trim(judan_csv), status="unknown")
			ss = 0.
			write(10,100) "#id", ",", "kp", ",", "xx", ",", "yy", ",", "ss", "," &
				, (trim(cname(n))//"_InitialBed(m)" , "," &
				,  trim(cname(n))//"_Elevation(m)" , "," &
				,  trim(cname(n))//"_WaterSurface(m)", "," &
				,  trim(cname(n))//"_VelocityX(m/s)", "," &
				,  trim(cname(n))//"_VelocityY(m/s)", "," &
				, n=0, nmax-1)
			
			id =1
			write(10,101) id, ",", judan(id)%name, ",", judan(id)%xx, ",", judan(id)%yy, "," , slen - ss, "," &
								, (zz0(n, judan(id)%np), "," &
								,  elv(n, judan(id)%np), "," &
								,  srf(n, judan(id)%np), "," &
								,  vx(n, judan(id)%np), "," &
								,  vy(n, judan(id)%np), "," &
								, n = 0, nmax-1)
			do id = 2, mmax
	
				ds = sqrt((judan(id)%xx - judan(id-1)%xx)**2. + (judan(id)%yy - judan(id-1)%yy)**2.)
				ss = ss + ds
				write(10,101) id, ",", judan(id)%name, ",", judan(id)%xx, ",", judan(id)%yy, ",", slen - ss, "," &
								, (zz0(n, judan(id)%np), "," &
								,  elv(n, judan(id)%np), "," &
								,  srf(n, judan(id)%np), "," &
								,  vx(n, judan(id)%np), "," &
								,  vy(n, judan(id)%np), "," &
								, n = 0, nmax-1)
		
			end do
			close(10)
			!----------------------------------------------------------------------
		
	100 format(a, a, a, 100(a, a))
	101 format(i5, a, a, 100(a, f20.8))		
	102 format(i5, 100(a, f20.8))		
		
			!----------------------------------------------------------------------
			!���f�}�p�f�[�^�̏o��
			!oudan_csv = 'kp=+0000.0_t=00000000.00(sec).csv'
			write(oudan_csv(14:24),'(f10.2)') tt
			do ii=14,24
				if(oudan_csv(ii:ii) == ' ') oudan_csv(ii:ii)='0'
			end do
		
			do id = 1, mmax
				read(odan(id)%name,*) skp
				write(oudan_csv(5:10),'(f6.1)') abs(skp)
				do ii=5,10
					if(oudan_csv(ii:ii) == ' ') oudan_csv(ii:ii)='0'
				end do
				if(skp<0)then
					oudan_csv(4:4) = "-"
				else
					oudan_csv(4:4) = "+"
				end if
			
				open(10, file= trim(juoudan_folder)//"\oudan\"//trim(oudan_csv), status="unknown")
				write(10,"(i10, a, a)") id, ",", "kp"//trim(odan(id)%name)
			
				ss = 0.
				jp =1
				np = odan(id)%np(jp)
				write(10,100) "#id", ",", "xx", ",", "yy", ", ", "ss", "," &
					, (trim(cname(n))//"_InitialBed(m)", "," &
					,  trim(cname(n))//"_Elevation(m)", "," &
					,  trim(cname(n))//"_WaterSurface(m)", "," &
					,  trim(cname(n))//"_VelocityX(m/s)", "," &
					,  trim(cname(n))//"_VelocityY(m/s)", "," &
					, n=0, nmax-1)
				
				write(10, 102) jp, ",", xx(np), ",", yy(np), ",", ss, "," &
					, (zz0(n,np), "," &
					,  elv(n,np), "," &
					,  srf(n,np), "," &
					,  vx(n,np), "," &
					,  vy(n,np), "," &
					, n=0, nmax-1)
			
				do jp = 2, odan(id)%jpmax

					np0 = odan(id)%np(jp-1)
					np = odan(id)%np(jp)

					ds = sqrt((xx(np) - xx(np0))**2. + (yy(np) - yy(np0))**2.)
					ss = ss + ds
					write(10, 102) jp, ",", xx(np), ",", yy(np), ",", ss, "," &
						, (zz0(n,np), "," &
						,  elv(n,np), "," &
						,  srf(n,np), "," &
						,  vx(n,np), "," &
						,  vy(n,np), "," &
						, n=0, nmax-1)
			
				end do
				close(10)
				!----------------------------------------------------------------------
			
			end do
		
		end if

   end do

	return
	
999	write(*,*) "CAT2D can be applied only to the same grid size condition."
		stop
		
998	write(*,*) "Unstructured is specified as Grid type condition."
		write(*,*) "But following CGNS file does not have Unstructured Grid."
		write(*,*) "Please check your calculation condition."
		stop

997	write(*,*) "CAT2D can be applied only to the same number of output condition."
		stop

		
end subroutine