!**********************************************************************
!
! GUI�Ŏw�肵��parametername.csv ����
!	��r�p�����[�^�����擾���A�O���[�o���ϐ��Ƃ��Đݒ肷��
!
!**********************************************************************
   subroutine set_comparison_parametername
   use common_params
   implicit none
   integer :: i, ii, n, icount, l
   character(len=strmax) :: w, a
!----------------------------------------------------------------------
   
	pname = ""
	
   open(1, file=trim(pname0), status="old", err=999 )  
   read (1, '()')       ! �w�b�_�s�̓ǂݔ�΂�
   
   icount = 0
   do
      read (1, *, end=100) ii, w  ! �t�@�C���I�[�Ȃ��999�ɔ��
      icount = icount + 1
   end do
100 continue   
   
   !��r�p�����[�^�̐�
   num_p = icount
   allocate(pname(0:nmax-1, 1:num_p))
   
   rewind (1)           !�@�t�@�C���̓��ɖ߂�
   read (1, '()')       !�@�w�b�_�s�̓ǂݔ�΂�
   do i = 1, num_p
      read(1, *, err = 998) ii, (pname(n,i), n=0, nmax-1)
   end do
	
	do i= 1, num_p
		do n= 0, nmax-1
			l = len_trim(pname(n,i))
			!l = len(a)
			if(l == strmax )then
				write(*,*) "Error : Blank is used for a parameter name."
				write(*,*) "Please check your specified file."
				stop
			end if
		end do
	end do
   
   close(1)
	return
   
999 continue
	write(*,*) "Error : File path is wrong. Please check it."
	write(*,*) "File Path = ", trim(pname0)
	stop
	
998 continue
	write(*,*) "Error : Parameter File format is invalid. Please check it."
	stop
	
end subroutine