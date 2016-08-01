subroutine structured_grid_analysis
	use common_params
	implicit none
   include 'cgnslib_f.h'
	
   !カウンタ
	integer:: n				!ケース数ループカウンタ
	integer:: m, id		!横断面ループカウンタ
	integer:: jp			!横断面・測点ループカウンタ
	integer:: i, j			!格子点ループカウンタ
	integer:: l				!出力数ループカウンタ
	integer:: ii			!出力ファイル名作成用
	!integer:: i, j, k, l, m, n
	
	!最短格子点の探索
	double precision :: ds, dsmin, ss
	
	!キロポスト値
	double precision :: skp
   
	!エラー判定値
	integer :: ier
	
	!各ケースの計算結果出力回数、lmax0はBaseCase
	integer :: lmax0, lmax
	
	!比較格子サイズ・格子
   integer:: isize0, jsize0, isize, jsize
   double precision, dimension(:,:), allocatable:: xx, yy
	
   !時間
   double precision :: tt
	
	!河床変動計算フラグ
	integer :: ibedcal
    
   !格子物理変数
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

    
   !演算値
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
	
	!出力ファイル名
	character(len=22) :: judan_csv
	character(len=33) :: oudan_csv
	
	!work用
	double precision, dimension(:), allocatable:: wtmp
	integer :: iip, jjp, iip0, jjp0
!================================================================================
	
	!--------------------------------------------------
	!出力用フォルダ、ファイル名の設定
	if(irsd == 1)then
		!フォルダ有無確認
		call system("mkdir "//trim(juoudan_folder)//"\judan")
		call system("del /Q "//trim(juoudan_folder)//"\judan\*.csv")
		call system("mkdir "//trim(juoudan_folder)//"\oudan")
		call system("del /Q "//trim(juoudan_folder)//"\oudan\*.csv")
	
		!出力ファイル名デフォルト
		judan_csv = 't=00000000.00(sec).csv'
		oudan_csv = 'kp=+0000.0_t=00000000.00(sec).csv'
	end if
	
	!----------------------------------------------------------------------
   !BaseCaseの計算結果から
	! 計算結果出力回数：lmaxの設定  →　同じ出力回数でないと比較しない　ただし、同時刻か否かは判定していない
   ! 計算結果出力回数：lmaxの設定  →　少ないほうに合わせる。ただし、同時刻か否かは判定していない
	
   call cg_iric_read_sol_count_mul_f(idf(0), lmax0, ier)
   
	!出力回数の異なる比較ケースがないかエラーチェック
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
   !BaseCaseの計算結果から
	! 比較の基礎となる格子:isize, jsize, xx, yyを設定する　→　cgnsファイルに出力
   call cg_iric_gotogridcoord2d_mul_f(idf(0), isize0, jsize0, ier)
	if(ier /= 0) then
		write(*,*) "***** Error : Reading File *****"
		write(*,*) "File = ", trim(fname(0))
		goto 998
	end if
	
	!格子数の異なる比較ケースがないかエラーチェック
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

   !比較する計算対象格子を確定
   isize = isize0
   jsize = jsize0
	write(*,*) "Grid isize", isize
	write(*,*) "Grid jsize", jsize
	
   allocate(xx(1:isize, 1:jsize), yy(1:isize, 1:jsize))
   call cg_iric_getgridcoord2d_mul_f(idf(0), xx, yy, ier)
   call cg_iric_writegridcoord2d_mul_f(idf0, isize, jsize, xx, yy, ier)
	
	
	!----------------------------------------------------------------------
	!縦断図用　
	!→　中心点、および、各横断測点に最も近い格子点idを設定
	if(irsd == 1)then
	
		!河道中心点に最も近い格子点番号
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
	!横断図用　
	!→　中心点、および、各横断測点に最も近い格子点idを設定
	if(irsd == 1)then
	
		!各横断測点に最も近い格子点番号
		do m = 1, mmax  !断面ループ
			do jp = 1, odan(m)%jpmax  !測点ループ
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
   
	!差分値の出力
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