!**********************************************************************
!
!	�e�P�[�X����v�Z���ʂ�Ǎ��������l�쐬���o��
!
!**********************************************************************
subroutine rescomp(lmax, isize, jsize, xx, yy)
	use common_params
	implicit none
	include 'cgnslib_f.h'
	include 'iriclib_f.h'
    
	integer, intent(in) :: lmax, isize, jsize
	double precision, dimension(1:isize,1:jsize), intent(in):: xx, yy
	
	integer, dimension(:,:,:,:), allocatable :: ival, delt_ival
	double precision, dimension(:,:,:,:), allocatable :: rval, delt_rval
	
	!�i�q�����ϐ�
	double precision, dimension(:,:,:), allocatable::  zz0
   double precision, dimension(:,:,:), allocatable::  sum_v, delt_sum_v
	
	!�o�̓t�@�C����
	character(len=22) :: judan_csv
	character(len=33) :: oudan_csv
    
    !�p�����[�^���̕������`�F�b�N�p
    character(len=strmax) :: wbuf
	
	integer :: i,j,k,l,m,n, ii, iip0, jjp0, iip, jjp
	integer :: ier, id, jp
	double precision :: tt, ds, ss, skp
!======================================================================
	
	!----------------------------------------------------------------------
	allocate(ival(0:nmax-1, 1:num_p, 1:isize, 1:jsize))
	allocate(rval(0:nmax-1, 1:num_p, 1:isize, 1:jsize))
	allocate(delt_ival(0:nmax-1, 1:num_p, 1:isize, 1:jsize))
	allocate(delt_rval(0:nmax-1, 1:num_p, 1:isize, 1:jsize))

	allocate(zz0(0:nmax-1, 1:isize, 1:jsize))
   allocate(sum_v(0:nmax-1, 1:isize, 1:jsize), delt_sum_v(0:nmax-1, 1:isize, 1:jsize))
	
	!----------------------------------------------------------------------
	!�o�͗p�t�H���_�A�t�@�C�����̐ݒ�
	if(juoudan_output == 0 .and. irsd == 1)then
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
	!�v�Z���ʓǍ��J�n
   do l = 1, lmax
	 
		!----------------------------------------------------------------------
      !BaseCase���玞�Ԃ�ǂݍ���
		call cg_iric_read_sol_time_mul_f(idf(0), l, tt, ier)		
		
		!----------------------------------------------------------------------
		!�v�Z���ʂ̓Ǎ�
      do n = 0, nmax-1
			do m = 1, num_p
				call cg_iric_read_sol_real_mul_f(idf(n), l, pname(n,m), rval(n,m,:,:), ier)
			end do
		end do
		
		!----------------------------------------------------------------------
		!�����͏���ۑ�
		if(l==1)then
			do n = 0, nmax-1
				do i = 1, isize
					do j = 1, jsize
						zz0(n,i,j) = rval(n,1,i,j)
					end do
				end do
			end do
		end if
      
      !----------------------------------------------------------------------
		!�����x�N�g���̐�Βl
		do n = 0, nmax-1
			do i = 1, isize
				do j = 1, jsize
					sum_v(n,i,j) = sqrt(rval(n,3,i,j)**2. +rval(n,4,i,j)**2.)
				end do
			end do
		end do
		
      !�����x�N�g���̍���
      do n = 1, nmax-1
			do i = 1, isize
				do j = 1, jsize
					delt_sum_v(n,i,j) = sum_v(n,i,j) - sum_v(0,i,j)
				end do
			end do
      end do
    
    ! ���[�U��GUI��� "STOP" �{�^���������Ď��s���L�����Z���������m�F
    call iric_check_cancel_f(ier)
    if (ier == 1) then
        write(*,*) "Solver is stopped because the STOP button was clicked."
        stop
    end if
      
    !
    ! gui��cgns�t�@�C����Ǎ������ۂ��𔻒�
    !
    do
        call iric_check_lock_f(fname0, ier)
        if (ier == 1) then
            call sleep(1)
        elseif (ier == 0) then  !�Ǎ����łȂ����do���[�v�𔲂���
            exit
        end if
    end do
    call iric_write_sol_start_f(fname0, ier)
      
	!----------------------------------------------------------------------
      !�����l�̎Z�o�l
		do m = 1, num_p
			do n = 1, nmax-1
				do j=1, jsize
					do i=1, isize
						delt_rval(n,m,i,j) = rval(n,m,i,j) - rval(0,m,i,j)
					end do
				end do
			end do
		end do
		
		!----------------------------------------------------------------------
		!���Ԃ��o��
		call cg_iric_write_sol_time_mul_f(idf0, tt, ier)
      write(*,*) "time = ", tt, "[sec]"
		
		!----------------------------------------------------------------------
		!��r�ł���悤�@���̂܂܂��ׂĂ��o��
		do n = 0, nmax-1
			do m = 1, num_p
                wbuf=""
                wbuf = trim(cname(n))//"_"//trim(pname(n,m))
                !write(*,*) len_trim(wbuf)
                if(len_trim(wbuf)>30) then 
                    wbuf = wbuf(1:30)
                    
                    !write(*,"(a)") "Warning:Output parameter's name is too long for the cgns format."
                    !write(*,"(a)") "So I would like to reccomend to change the case's name into a short one."
                    !call iric_write_sol_end_f(fname0, ier)
                    !stop
                end if
				call cg_iric_write_sol_real_mul_f(idf0, wbuf(1:30), rval(n,m,:,:), ier)
			end do
      end do
		
		!----------------------------------------------------------------------
		!�����l���o��
		do n = 1, nmax-1
			do m = 1, num_p
                wbuf=""
                !wbuf = trim(cname(n))//"-"//trim(cname(0))//"_"//trim(pname(n,m))//"_"
                wbuf = trim(cname(n))//"-"//trim(cname(0))//"_"//trim(pname(n,m))
                if(len_trim(wbuf)>30) then 
                    
                    wbuf = wbuf(1:30)
                    
                    !write(*,"(a)") "Warning:Output parameter's name is too long for the cgns format."
                    !write(*,"(a)") "So I would like to reccomend to change the case's name into a short one."
                    !call iric_write_sol_end_f(fname0, ier)
                    !stop
                end if
                call cg_iric_write_sol_real_mul_f(idf0, wbuf(1:30), delt_rval(n,m,:,:), ier)
            end do
         
         !���������̍����l
            wbuf=""
            wbuf = trim(cname(n))//"-"//trim(cname(0))//"_SumVelocity"
                if(len_trim(wbuf)>30) then 
                    
                    wbuf = wbuf(1:30)
                    
                    !write(*,"(a)") "Warning:Output parameter's name is too long for the cgns format."
                    !write(*,"(a)") "So I would like to reccomend to change the case's name into a short one."
                    
                end if
		   call cg_iric_write_sol_real_mul_f(idf0, wbuf(1:30), delt_sum_v(n,:,:), ier)
        end do
      
        !HDD�֓f���o��
        call cg_iric_flush_f(fname0, idf0, ier)
      
        !�o�͏����I��
        call iric_write_sol_end_f(fname0, ier)
		
		!----------------------------------------------------------------------
		!�c�f�f�[�^�̏o�́@���@csv�`���ŏo��(judan_csv = 't=00000000.00(sec).csv')
		if(juoudan_output == 0 .and. irsd == 1)then

			write(judan_csv(3:13),'(f10.2)') tt
			do ii=3,13
				if(judan_csv(ii:ii) == ' ') judan_csv(ii:ii)='0'
			end do
		
			open(10, file=trim(juoudan_folder)//"\judan\"//judan_csv, status="unknown")
			
			!Header
			write(10,100) "#id", ",", "kp", ",", "xx", ",", "yy", ",", "ss", "," &
				, (trim(cname(n))//"_InitialBed(m)" , "," &
				,  trim(cname(n))//"_Elevation(m)" , "," &
				,  trim(cname(n))//"_WaterSurfaceElevation(m)", "," &
				,  trim(cname(n))//"_VelocityX(m/s)", "," &
				,  trim(cname(n))//"_VelocityY(m/s)", "," &
				, n=0, nmax-1)
	
			!�o��
			ss = 0
			do id = 1, mmax
				if(id>1)then
					ds = sqrt((judan(id)%xx - judan(id-1)%xx)**2. + (judan(id)%yy - judan(id-1)%yy)**2.)
					ss = ss + ds	
				end if
				
				write(10,101) id, ",", trim(judan(id)%name), ",", judan(id)%xx, ",", judan(id)%yy, ",", slen - ss, "," &
								, (zz0(   n, judan(id)%ip, judan(id)%jp), "," &
								,  rval(   n, 1, judan(id)%ip, judan(id)%jp), "," &
								,  rval(   n, 2, judan(id)%ip, judan(id)%jp), "," &
								,  rval(   n, 3, judan(id)%ip, judan(id)%jp), "," &
								,  rval(   n, 4, judan(id)%ip, judan(id)%jp), "," &
								, n = 0, nmax-1)				
				
			end do
			close(10)
		end if
		
		
		!----------------------------------------------------------------------
		!���f�f�[�^�̏o�� ���@csv�`���ŏo��(!oudan_csv = 'kp=+0000.0_t=00000000.00(sec).csv')
		!----------------------------------------------------------------------
		if(juoudan_output == 0 .and. irsd == 1)then
			
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
			
				!�o�̓t�@�C�����J��
				open(10, file= trim(juoudan_folder)//"\oudan\"//trim(oudan_csv), status="unknown")
				!write(10,"(i10, a, a)") id, ",", "kp"//trim(odan(id)%name)
				
				!Header
				write(10,100) "#id", ",", "xx", ",", "yy", ", ", "ss", "," &
					, (trim(cname(n))//"_InitialBed(m)", "," &
					,  trim(cname(n))//"_Elevation(m)", "," &
					,  trim(cname(n))//"_WaterSurface(m)", "," &
					,  trim(cname(n))//"_VelocityX(m/s)", "," &
					,  trim(cname(n))//"_VelocityY(m/s)", "," &
					, n=0, nmax-1)
					
				
				!�o��
				ss = 0.
				do jp = 1, odan(id)%jpmax
					iip  = odan(id)%ip(jp)
					jjp  = odan(id)%jp(jp)
					
					if(jp > 1)then
						iip0 = odan(id)%ip(jp-1)
						jjp0 = odan(id)%jp(jp-1)
						ds = sqrt((xx(iip, jjp) - xx(iip0, jjp0))**2. + (yy(iip, jjp) - yy(iip0, jjp0))**2.)
						ss = ss + ds
					end if
					write(10, 102) jp, ",", xx(iip, jjp), ",", yy(iip, jjp), ",", ss, "," &
						, (zz0(n, iip, jjp), "," &
						,  rval(n, 1, iip, jjp), "," &
						,  rval(n, 2, iip, jjp), "," &
						,  rval(n, 3, iip, jjp), "," &
						,  rval(n, 4, iip, jjp), "," &
						, n=0, nmax-1)
			
				end do

				close(10)
			end do
		end if
		
	end do
	
	return
	
	100 format(a, a, a, 100(a, a))			!�c�f�}�w�b�_�[�p
	101 format(i5, a, a, 100(a, f20.8))		!���f�}�w�b�_�[�p
	102 format(i5, 100(a, f20.8))				!�f�[�^�o�͗p
	
end subroutine