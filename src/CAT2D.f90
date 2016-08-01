program cat2d
	use common_params
   implicit none
   include 'cgnslib_f.h'
	include 'iriclib_f.h'

	integer :: ibase, izone
	integer :: zonetype0, zonetype
	
	!�͐쑪�ʃf�[�^�p
	integer :: num_geo
	integer :: geoid
	integer :: geotype
	
	!�J�E���^
   integer:: n
	
   !�G���[�t���O
   integer:: ier
   
   !====================================================================================================

	!--------------------------------------------------
   ! get cgns file name
   if (nargs() == 2) then
      call getarg(1, fname0, ier)
   else
      write(*,*) "Input File not specified."
      stop
   end if
	
   !--------------------------------------------------
   ! Open CGNS file
   call cg_open_f(trim(fname0), CG_MODE_MODIFY, idf0, ier)
   if (ier /=0) STOP "*** Open error of CGNS file ***"
   call cg_iric_init_f(idf0, ier)
	
   !--------------------------------------------------
   !��r�P�[�X���̎擾�@���@�e�P�[�X�̕ϐ����A���P�[�g 
	! ���@�x�[�X�P�[�X�{�R�P�[�X�܂Ŕ�r�ł���悤definition.xml�Ő��䂵�Ă���i2016.03.02�j
	! 0:BaseCase
	! 1:Case1
	! 2:Case2
	! 3:Case3
   call cg_iric_read_integer_mul_f(idf0, "nmax", nmax, ier)
	allocate(idf(0:nmax-1), fname(0:nmax-1), cname(0:nmax-1))
	
	!--------------------------------------------------
   !�t�@�C�����̎擾
	!�K�{
   call cg_iric_read_string_mul_f(idf0, "cname0", cname(0), ier)
   call cg_iric_read_string_mul_f(idf0, "fname0", fname(0), ier)
   call cg_iric_read_string_mul_f(idf0, "cname1", cname(1), ier)
   call cg_iric_read_string_mul_f(idf0, "fname1", fname(1), ier)
   
	!�C��
   if(nmax>2) call cg_iric_read_string_mul_f(idf0, "cname2", cname(2), ier)
   if(nmax>2) call cg_iric_read_string_mul_f(idf0, "fname2", fname(2), ier)
   if(nmax>3) call cg_iric_read_string_mul_f(idf0, "cname3", cname(3), ier)
   if(nmax>3) call cg_iric_read_string_mul_f(idf0, "fname3", fname(3), ier)
   if(nmax>4) call cg_iric_read_string_mul_f(idf0, "cname4", cname(4), ier)
   if(nmax>4) call cg_iric_read_string_mul_f(idf0, "fname4", fname(4), ier)
   if(nmax>5) call cg_iric_read_string_mul_f(idf0, "cname5", cname(5), ier)
   if(nmax>5) call cg_iric_read_string_mul_f(idf0, "fname5", fname(5), ier)
   
   !--------------------------------------------------
   !��r�p�����[�^�����擾����
   call cg_iric_read_string_mul_f(idf0, "pname0", pname0, ier)
   
   !��r�p�����[�^�[����ݒ肷��
   call set_comparison_parametername
   
   !--------------------------------------------------
   !��r����t�@�C�����J��
   do n = 0, nmax-1
      call cg_open_f(trim(fname(n)), CG_MODE_READ, idf(n), ier)
      if(ier /= 0)then
			write(*,"(a)") "***** File open Error. *****"
			write(*,"(a)") "CGNS file Path is missing. Please check your calculation condition."
			write(*,"(a)") "Specified CGNS file Path is ", trim(fname(n))
			stop
		end if
		
		call cg_iric_initread_f(idf(n), ier)
	end do
	
   !--------------------------------------------------
   !�c�E���f�f�[�^�o�̓t���O
   call cg_iric_read_integer_mul_f(idf0, "juoudan_output", juoudan_output, ier)
   
   !�c�E���f�f�[�^�o�̓t�H���_
   call cg_iric_read_string_mul_f(idf0, "juoudan_folder", juoudan_folder, ier)
   
	!--------------------------------------------------
   !�͐쑪�ʃf�[�^������΂����Ǎ���
   irsd = -1
   if(juoudan_output == 0)then
	   call cg_iric_read_geo_count_mul_f(idf0, "Elevation", num_geo, ier)
		
		if(num_geo > 2)then
			write(*,"(a)") "***** Warning *****"
			write(*,"(a)") "You've imported more than two River Survey Data."
			write(*,"(a)") "Upper position of River Survey Data in the Object Browser is effective"
			write(*,"(a)") "to output Longitudinal and Cross Sectional."
			write(*,"(a)") "***** Warning *****"
		end if
		
	   do geoid = 1, num_geo
		   call cg_iric_read_geo_filename_mul_f(idf0, "Elevation", geoid, gname, geotype, ier)	
		   if(geotype == iRIC_GEO_RIVERSURVEY)then
			   irsd = 1
			   exit
		   end if
	   end do
	   
	   if(irsd == 1)then
		   !write(*,"(a)") "RiverSurveyData = Yes"
		   call read_riversurvey
		
	   else
		   write(*,"(a)") "***** Error *****"
			write(*,"(a)") "You specified 'Valid' for the Longitudinal and Cross Sectional Data Output."
			write(*,"(a)") "But the RiverSurveyData is not imported."
			write(*,"(a)") "Please check your calculation condition."
			write(*,"(a)") "***** Error *****"
			stop
			
	   end if
   end if
	
	!--------------------------------------------------
	!�\���o�[�����m�F����
	call cg_gopath_f(idf(0), '/iRIC/SolverInformation', ier)
	!call cg_narrays_f(narrays, ier)
	call cg_array_read_f(1, solname0, ier)
	
	!--------------------------------------------------
	!Nays2DH��Mflow�ɂ̂ݑΉ����Ă���
	!if( &
	!	solname0(1:7) == "Nays2DH" .or. &
	!	solname0(1:8) == "Mflow_02" )then
	!	
	!	write(*,"(a, a)") "Solver Name = ", trim(solname0)
	!	
	!else
	!	write(*,"(a)") "***** Error *****"
	!	write(*,"(a)") "Current CAT2D is implemented for Nays2DH or Mflow_02."
	!	stop
	!end if
	
	!--------------------------------------------------
	!!�قȂ�\���o�[�͔�r�ł��Ȃ�
	!do n = 1, nmax-1
	!	call cg_gopath_f(idf(n), '/iRIC/SolverInformation', ier)
	!	call cg_array_read_f(1, solname, ier)
	!	
	!	if(solname /= solname0)then
	!		write(*,"(a)") "Current version of this tool cannot compare the different solver's result."
	!		stop
	!	end if
	!end do

	!--------------------------------------------------
	!�قȂ�i�q�\���̈قȂ�\���o�[�͔�r�ł��Ȃ�
	!�\���i�q or ��\���i�q�t���O���擾�@�\���i�q��zonetype=2�A��\���i�q��zonetype=3
	ibase = 1
	izone =1
	call cg_zone_type_f(idf(0), ibase, izone, zonetype0, ier)
	do n = 1, nmax-1
      call cg_zone_type_f(idf(0), ibase, izone, zonetype, ier)
		if(zonetype /= zonetype0)then
			write(*,*) "Different type of grid system can not be compared."
			stop
		end if
	end do

	if(zonetype0 == 2) then
		!�\���i�q��
		write(*,"(a)") "Grid Type = Structured Grid "
		call structured_grid_analysis
	
	elseif(zonetype0 == 3)then
		!��\���i�q��
		write(*,"(a)") "Grid Type = Unstructured Grid "
		call unstructured_grid_analysis
	
	else
		write(*,*) "Unknown grid type was specified. Please check your grid type in your calculation conditions."
		stop
	end if

   !--------------------------------------------------
   ! CGNS �t�@�C���̃N���[�Y
   call cg_close_f(idf0, ier)
   do n = 0, nmax-1
      call cg_close_f(idf(n), ier)
   end do
   
   write(*,*) "CGN Analysis tool normal end. Congraturation!!!"
   stop
	
end program