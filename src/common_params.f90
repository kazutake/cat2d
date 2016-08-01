module common_params
	
	integer, parameter :: strmax = 1024
	
	type PointSet
		character(len=strmax) :: name
		double precision :: xx, yy
		integer :: ip, jp		!構造格子用
		integer :: np			!非構造格子用
	end type
	
	type CrossSection
		character(len=strmax) :: name
		integer :: jpmax
		double precision, dimension(:), allocatable :: xx, yy
		integer, dimension(:), allocatable :: ip, jp		!構造格子用
		integer, dimension(:), allocatable :: np			!非構造格子用
	end type
	
	!比較ケース数
   integer :: nmax
   
   !縦横断図出力フラグ、フォルダ名
   integer :: juoudan_output
   character(len=strmax) :: juoudan_folder
   
	!ソルバー名
	character(len=strmax) :: solname0, solname
	
   !ファイル名
   integer :: idf0
   character(len=strmax) :: fname0
   integer, dimension(:), allocatable :: idf
   character(len=strmax), dimension(:), allocatable :: cname, fname
   
   !変数名称
   integer :: num_p
   character(len=strmax) :: pname0
	character(len=strmax), dimension(:,:), allocatable :: pname
   
	!河川測量データ
	integer :: irsd							!フラグ
	character(len=1024) :: gname			!データファイル名
	integer :: mmax							!横断面数
	type(PointSet), dimension(:), allocatable :: judan
	type(CrossSection), dimension(:), allocatable :: odan
	double precision :: slen				!中心線に沿った流下方向距離
	
end module