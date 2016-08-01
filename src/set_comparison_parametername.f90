!**********************************************************************
!
! GUIで指定したparametername.csv から
!	比較パラメータ名を取得し、グローバル変数として設定する
!
!**********************************************************************
   subroutine set_comparison_parametername
   use common_params
   implicit none
   integer :: i, ii, n, icount, l
   character(len=strmax) :: w, a
	
	integer :: ipos1, ipos2,buf_len
	character(len=strmax) :: buf, row_id
!----------------------------------------------------------------------
   
	pname = ""
	
   open(1, file=trim(pname0), status="old", err=999 )  
   read (1, '()')       ! ヘッダ行の読み飛ばし
   
   icount = 0
   do
      read (1, *, err=100, end=100) ii, w  ! ファイル終端ならば999に飛ぶ
		icount = icount + 1
   end do
100 continue   
   
   !比較パラメータの数
   num_p = icount
   allocate(pname(0:nmax-1, 1:num_p))
   
   rewind (1)           !　ファイルの頭に戻る
   read (1, '()')       !　ヘッダ行の読み飛ばし
	
   do i = 1, num_p
      !read(1, *, err = 998) ii, (pname(n,i), n=0, nmax-1)
		
		read(1, "(a)", err = 998) buf
		
		buf_len = len_trim(buf)
		ipos1 = 1	
		n=-1
		do 
			buf_len = len_trim(buf)
			ipos2 = scan(buf, ",")
			if(n>=0)then
				if(ipos2==0) then
					ipos2 = buf_len
					if(buf(ipos1:ipos2-1)=="") exit
					pname(n,i) = buf(ipos1:ipos2)
				else
					if(buf(ipos1:ipos2-1)=="") exit
					pname(n,i) = buf(ipos1:ipos2-1)
				end if
			end if
			
			buf = buf(ipos2+1:len(buf))
			
			n = n+1
		end do 
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