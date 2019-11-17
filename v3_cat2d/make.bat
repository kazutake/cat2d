ifort 	    ..\libs\cgnsdll_x64_ifort.lib ^
	    ..\libs\iriclib_x64_ifort.lib ^
	    ..\src\common_params.f90 ^
	    ..\src\unst_grid.f90 ^
	    ..\src\for_structured.f90 ^
	    ..\src\for_unstructured.f90 ^
	    ..\src\length.f90 ^
	    ..\src\outcsv_judan.f90 ^
	    ..\src\read_riversurvey.f90 ^
	    ..\src\rescomp.f90 ^
	    ..\src\set_comparison_parametername.f90 ^
   	    ..\src\CAT2D.f90 ^
	     /O2 /MD /I ..\include -o cat2d.exe


rm *.obj
rm *.mod



