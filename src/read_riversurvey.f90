subroutine read_riversurvey
	use common_params
	implicit none
	integer :: ier
	
	integer :: rsid, id
	integer :: m
	
	!work用変数
	double precision :: xc, yc, xdir, ydir
	double precision, dimension(:), allocatable :: altpos, altheight
	integer, dimension(:), allocatable :: altactive
	double precision :: sl, coss, sinn
	double precision :: ss, ds
!======================================================================
	
	!ファイルを開く
	call iric_geo_riversurvey_open_f(gname, rsid, ier)
	
	!横断面数を取得
	call iric_geo_riversurvey_read_count_f(rsid, mmax, ier) 
	allocate(judan(1:mmax))
	allocate(odan( 1:mmax))
	
	!横断面数のループ
	do id = 1, mmax	
	
		!縦断図となるxy座標値＝横断線の中心点座標を取得
		call iric_geo_riversurvey_read_position_f(rsid, id, judan(id)%xx, judan(id)%yy, ier)
		call iric_geo_riversurvey_read_name_f(rsid, id, judan(id)%name, ier)
		
		!横断データの名前と測点数を取得
		call iric_geo_riversurvey_read_name_f(rsid, id, odan(id)%name, ier)
		call iric_geo_riversurvey_read_altitudecount_f(rsid, id, odan(id)%jpmax, ier)
		
		!work用変数のメモリ確保
		allocate(altpos(odan(id)%jpmax))
		allocate(altheight(odan(id)%jpmax))
		allocate(altactive(odan(id)%jpmax)) 
		
		!各測点の値をwork用変数に格納
		call iric_geo_riversurvey_read_position_f( rsid, id, xc, yc, ier)
		call iric_geo_riversurvey_read_direction_f(rsid, id, xdir, ydir, ier) 
		call iric_geo_riversurvey_read_altitudes_f(rsid, id, altpos, altheight, altactive, ier) 
			
		!グローバル変数へx,yの値として代入する
		allocate(odan(id)%xx(odan(id)%jpmax), odan(id)%yy(odan(id)%jpmax)) 
		allocate(odan(id)%ip(odan(id)%jpmax), odan(id)%jp(odan(id)%jpmax))
		allocate(odan(id)%np(odan(id)%jpmax))
		sl = sqrt(xdir**2. + ydir**2.)
		coss = xdir / sl
		sinn = ydir / sl
		do m = 1, odan(id)%jpmax
			odan(id)%xx(m) = xc + altpos(m) * coss
			odan(id)%yy(m) = yc + altpos(m) * sinn
		end do
		
		!work用変数のメモリ解放
		deallocate(altpos)
		deallocate(altheight)
		deallocate(altactive)
		
	end do
	
	!流路長を算定
	ss = 0.
	do id = 2, mmax
		ds = sqrt((judan(id)%xx - judan(id-1)%xx)**2. + (judan(id)%yy - judan(id-1)%yy)**2.)
		ss = ss + ds
	end do
	slen = ss
	
	return
	
	
100 format(a5, 3(a, a20))
101 format(i5, 3(a1, f20.8))

end subroutine

