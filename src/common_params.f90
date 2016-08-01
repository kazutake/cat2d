module common_params
	
	integer, parameter :: strmax = 1024
	
	type PointSet
		character(len=strmax) :: name
		double precision :: xx, yy
		integer :: ip, jp		!�\���i�q�p
		integer :: np			!��\���i�q�p
	end type
	
	type CrossSection
		character(len=strmax) :: name
		integer :: jpmax
		double precision, dimension(:), allocatable :: xx, yy
		integer, dimension(:), allocatable :: ip, jp		!�\���i�q�p
		integer, dimension(:), allocatable :: np			!��\���i�q�p
	end type
	
	!��r�P�[�X��
   integer :: nmax
   
   !�c���f�}�o�̓t���O�A�t�H���_��
   integer :: juoudan_output
   character(len=strmax) :: juoudan_folder
   
	!�\���o�[��
	character(len=strmax) :: solname0, solname
	
   !�t�@�C����
   integer :: idf0
   character(len=strmax) :: fname0
   integer, dimension(:), allocatable :: idf
   character(len=strmax), dimension(:), allocatable :: cname, fname
   
   !�ϐ�����
   integer :: num_p
   character(len=strmax) :: pname0
	character(len=strmax), dimension(:,:), allocatable :: pname
   
	!�͐쑪�ʃf�[�^
	integer :: irsd							!�t���O
	character(len=1024) :: gname			!�f�[�^�t�@�C����
	integer :: mmax							!���f�ʐ�
	type(PointSet), dimension(:), allocatable :: judan
	type(CrossSection), dimension(:), allocatable :: odan
	double precision :: slen				!���S���ɉ�����������������
	
end module