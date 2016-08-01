!====================================================================================================
!
!> �i�q�̓ǂݍ��݂�R�s�[���s���������s�����W���[��
!
!====================================================================================================

module grid
	use common_params

	implicit none

	!include 'cgnslib_f.h'
	!include 'iriclib_f.h'

	!> �i�q�_�̐�
	integer :: gridSize
	!> �e�x�W��
	double precision, dimension(:), allocatable :: roughness
	!> ci ���ǂݍ��܂�Ă����� 1, ���Ȃ������� 0
	integer :: ciIsSet
	!> Channel Index (HSI��͗p)
	double precision, dimension(:), allocatable :: ci

contains

	!------------------------------------------------------------------------------------------
	!> �i�q�̔z����R�s�[
	!------------------------------------------------------------------------------------------
	subroutine copyArrayUnderGC(cgnsIn, cgnsOut, baseId, zoneId, gridId, gridIdTo, arrayId)
		integer, intent(in) :: cgnsIn, cgnsOut
		integer, intent(in) :: baseId, zoneId, gridId, gridIdTo
		integer, intent(in) :: arrayId

		character(len = strmax) :: arrayName
		integer :: dataType
		integer :: numDims
		integer :: dims(3)
		integer :: dataSize

		double precision, dimension(:), allocatable :: arrayData

		integer :: i
		integer :: ier

		! �z��̃��^�����擾
		call cg_goto_f(cgnsIn, baseId, ier, 'Zone_t', zoneId, 'GridCoordinates_t', gridId, 'end')
		call cg_array_info_f(arrayId, arrayName, dataType, numDims, dims, ier)

		dataSize = 1
		do i = 1, numDims
			dataSize = dataSize * dims(i)
		end do

		! �z��̒l��ǂݍ���
		allocate(arrayData(dataSize))
		call cg_array_read_f(arrayId, arrayData, ier)

		! �z��̒l����������
		call cg_goto_f(cgnsOut, baseId, ier, 'Zone_t', zoneId, 'GridCoordinates_t', gridIdTo, 'end')
		call cg_array_write_f(arrayName, dataType, numDims, dims, arrayData, ier)
		deallocate(arrayData)

	end subroutine

	!------------------------------------------------------------------------------------------
	!> �i�q�̑S�z����R�s�[
	!------------------------------------------------------------------------------------------
	subroutine copyAllArraysUnderGC(cgnsIn, cgnsOut, baseId, zoneId, gridId, gridIdTo)
		integer, intent(in) :: cgnsIn, cgnsOut
		integer, intent(in) :: baseId, zoneId, gridId, gridIdTo

		integer :: arrayId, numArrays

		integer :: ier

		call cg_goto_f(cgnsIn, baseId, ier, 'Zone_t', zoneId, 'GridCoordinates_t', gridId, 'end')
		call cg_narrays_f(numArrays, ier)

		do arrayId = 1, numArrays
			call copyArrayUnderGC(cgnsIn, cgnsOut, baseId, zoneId, gridId, gridIdTo, arrayId)
		end do

	end subroutine

	!------------------------------------------------------------------------------------------
	!> Zone�� Elements_t �v�f���R�s�[
	!------------------------------------------------------------------------------------------
	subroutine copyElement(cgnsIn, cgnsOut, baseId, zoneId, elemId)
		integer, intent(in) :: cgnsIn, cgnsOut
		integer, intent(in) :: baseId, zoneId, elemId

		character(len = strmax) :: sectionName
		integer :: sectionType, start, end, nbndry, parent_flag, dataSize, parentData
		integer, dimension(:), allocatable :: data
		integer :: tmpElemId

		integer :: ier

		! section �̏���ǂݍ���
		call cg_section_read_f(cgnsIn, baseId, zoneId, elemId, sectionName, &
			sectionType, start, end, nbndry, parent_flag, ier)
		call cg_elementdatasize_f(cgnsIn, baseId, zoneId, elemId, dataSize, ier)
		allocate(data(dataSize))
		call cg_elements_read_f(cgnsIn, baseId, zoneId, elemId, data, parentData, ier)

		! section �̏�����������
		call cg_section_write_f(cgnsOut, baseId, zoneId, sectionName, sectionType, &
			start, end, nbndry, data, tmpElemId, ier)

	end subroutine

	!------------------------------------------------------------------------------------------
	!> Zone�̑S Elements_t �v�f���R�s�[
	!------------------------------------------------------------------------------------------
	subroutine copyAllElements(cgnsIn, cgnsOut, baseId, zoneId)
		integer, intent(in) :: cgnsIn, cgnsOut
		integer, intent(in) :: baseId, zoneId

		integer :: elemId, numSections

		integer :: ier

		call cg_nsections_f(cgnsIn, baseId, zoneId, numSections, ier)

		do elemId = 1, numSections
			call copyElement(cgnsIn, cgnsOut, baseId, zoneId, elemId)
		end do

	end subroutine

	!------------------------------------------------------------------------------------------
	!> �Y������ Zone ���Ȃ���΁A�R�s�[����B
	!------------------------------------------------------------------------------------------
	subroutine copyZoneIfNeeded(cgnsIn, cgnsOut, baseId, zoneId, zoneCopied)
		integer, intent(in) :: cgnsIn, cgnsOut
		integer, intent(in) :: baseId, zoneId
		integer, intent(out) :: zoneCopied

		character(len = strmax) :: zoneName
		integer :: zoneSize(9)
		integer :: zoneType
		integer :: tmpZoneId
		integer :: numZones

		integer :: ier

		call cg_nzones_f(cgnsOut, baseId, numZones, ier)
		if (numZones .ne. zoneId - 1) then
			! Zone ���R�s�[����ׂ����ł͂Ȃ�
			zoneCopied = 0
			return
		end if

		! �܂��A�R�s�[���� Zone �̏���ǂݍ���
		call cg_zone_read_f(cgnsIn, baseId, zoneId, zoneName, zoneSize, ier)
		call cg_zone_type_f(cgnsIn, baseId, zoneId, zoneType, ier)

		! ���������R�s�[��ɏ�������
		call cg_zone_write_f(cgnsOut, baseId, zoneName, zoneSize, zoneType, tmpZoneId, ier)
		zoneCopied = 1

	end subroutine

	!------------------------------------------------------------------------------------------
	!> �i�q���R�s�[����B
	!------------------------------------------------------------------------------------------
	subroutine copyGrid(cgnsIn, cgnsOut, baseId, zoneId, gridId, gridIdTo, gridCopied)
		integer, intent(in) :: cgnsIn, cgnsOut
		integer, intent(in) :: baseId, zoneId, gridId, gridIdTo
		integer, intent(out) :: gridCopied

		character(len = strmax) :: coordName
		integer :: numGrids
		integer :: tmpGridId

		integer :: ier

		call cg_ngrids_f(cgnsIn, baseId, zoneId, numGrids, ier)
		
		if (numGrids < gridId) then
			! �R�s�[����i�q���Ȃ�
			gridCopied = 0
			return
		end if

		call cg_ngrids_f(cgnsOut, baseId, zoneId, numGrids, ier)
		
		if (numGrids /= gridIdTo - 1) then
			! �i�q���R�s�[����ׂ����ł͂Ȃ�
			gridCopied = 0
			return
		end if

		! GridCoordinates_t �̏���ǂݍ���
		call cg_grid_read_f(cgnsIn, baseId, zoneId, gridId, coordName, ier)
		! GridCoordinates_t ���쐬
		call cg_grid_write_f(cgnsOut, baseId, zoneId, coordName, tmpGridId, ier)

		! GridCoordinate_t �ȉ��̑S�Ă� �z����R�s�[
		call copyAllArraysUnderGC(cgnsIn, cgnsOut, baseId, zoneId, gridId, gridIdTo)

		gridCopied = 1

	end subroutine

	!------------------------------------------------------------------------------------------
	!> ���O�����Ċi�q���R�s�[����B
	!------------------------------------------------------------------------------------------
	subroutine copyGridAs(cgnsIn, cgnsOut, baseId, zoneId, gridId, gridIdTo, newName, gridCopied)
		integer, intent(in) :: cgnsIn, cgnsOut
		integer, intent(in) :: baseId, zoneId, gridId, gridIdTo
		character(len = *), intent(in) :: newName
		integer, intent(out) :: gridCopied

		character(len = strmax) :: coordName
		integer :: numGrids
		integer :: tmpGridId

		integer :: ier

		call cg_ngrids_f(cgnsIn, baseId, zoneId, numGrids, ier)
		
		if (numGrids < gridId) then
			! �R�s�[����i�q���Ȃ�
			gridCopied = 0
			return
		end if

		call cg_ngrids_f(cgnsOut, baseId, zoneId, numGrids, ier)
		
		if (numGrids /= gridIdTo - 1) then
			! �i�q���R�s�[����ׂ����ł͂Ȃ�
			gridCopied = 0
			return
		end if

		! GridCoordinates_t �̏���ǂݍ���
		call cg_grid_read_f(cgnsIn, baseId, zoneId, gridId, coordName, ier)
		! GridCoordinates_t ���쐬
		call cg_grid_write_f(cgnsOut, baseId, zoneId, newName, tmpGridId, ier)

		! GridCoordinate_t �ȉ��̑S�Ă� �z����R�s�[
		call copyAllArraysUnderGC(cgnsIn, cgnsOut, baseId, zoneId, gridId, gridIdTo)

		gridCopied = 1

	end subroutine

end module

!------------------------------------------------------------------------------------------
!> �i�q��ǂݍ��݁A�i�q�_�̐��� gridSize �ɏo�͂���
!------------------------------------------------------------------------------------------
subroutine gridLoad(cgnsIn, cgnsOut)
	use common_params
	use grid

	implicit none

	integer, intent(in) :: cgnsIn, cgnsOut
	integer :: i, j, ii, jj, index
	integer :: firstdata
	integer :: ier
	integer :: baseId, zoneId
	character(len = strmax) :: zoneName
	integer :: zoneSize(9)
	integer :: zoneType
	double precision, dimension(:,:), allocatable :: roughness_cell !< �e�x�W�� (�Z��)
	double precision :: tmp_roughness

	! BASE, ZONE �͊�{�I��1

	baseId = 1
	zoneId = 1

	! �i�q����ǂݍ���
	call cg_zone_read_f(cgnsIn, baseId, zoneId, zoneName, zoneSize, ier)
	call cg_zone_type_f(cgnsIn, baseId, zoneId, zoneType, ier)

	! �i�q�_�̐��� gridSize �Ɋi�[
	if (zoneType .eq. 2) then
		gridSize = zoneSize(1) * zoneSize(2)
	else 
		gridSize = zoneSize(1)
	end if

	! Roughness ��ǂݍ���
	allocate(roughness(gridSize))
	roughness = 0

	if (zoneType .eq. 2) then
		! �܂��Aroughness_cell �̓ǂݍ��݂����݂�B
		allocate(roughness_cell(zoneSize(1) - 1, zoneSize(2) - 1))
		call cg_iric_read_grid_real_node_mul_f(cgnsIn, 'roughness_cell', roughness_cell, ier)
		if (ier .eq. 0) then
			! �ǂݍ��݂����������B�Z���Œ�`����Ă���̂ŁA���̊i�q�_���͂�4�̃Z���ł�
			! �e�x�W���̍ő�l���̗p����
			do i = 1, zoneSize(1)
				do j = 1, zoneSize(2)
					firstdata = 1
					do ii = -1, 0
						do jj = -1, 0
							tmp_roughness = -1
							if (i + ii >= 1 .and. i + ii <= zoneSize(1) - 1 .and. &
							    j + jj >= 1 .and. j + jj <= zoneSize(2) - 1) then
								tmp_roughness = roughness_cell(i + ii, j + jj)
							end if
							if (tmp_roughness /= -1) then
								index = i + (j - 1) * (zoneSize(1))
								if (firstdata == 1 .or. tmp_roughness > roughness(index)) then
									roughness(index) = tmp_roughness
								end if
								firstdata = 0
							end if
						end do
					end do
				end do
			end do
		else
			! �ǂݍ��݂����s�����B�i�q�_�Œ�`����Ă���Ɖ���
			call cg_iric_read_grid_real_node_mul_f(cgnsIn, 'Roughness', roughness, ier)
			if (ier .ne. 0) then
				! �ǂݍ��݂Ɏ��s������A "roughness" ������
				call cg_iric_read_grid_real_node_mul_f(cgnsIn, 'roughness', roughness, ier)
			end if
		end if
	else 
		! roughness �́A�i�q�_�Œ�`����Ă���Ɖ���
		call cg_iric_read_grid_real_node_mul_f(cgnsIn, 'Roughness', roughness, ier)
		if (ier .ne. 0) then
			! �ǂݍ��݂Ɏ��s������A "roughness" ������
			call cg_iric_read_grid_real_node_mul_f(cgnsIn, 'roughness', roughness, ier)
		end if
		if (ier .ne. 0) then
			! �ǂݍ��݂Ɏ��s������A "initial_cd" ������ (SToRM �p)
			call cg_iric_read_grid_real_node_mul_f(cgnsIn, 'initial_cd', roughness, ier)
		end if
	end if

	! CI �̒l��ǂݍ���
	call loadCI(cgnsOut, ciIsSet, ci)

	! �Ō�ɁA�K�v������Ίi�q�� cgnsIn ���� cgnsOut �ɃR�s�[����
	call gridCopyIfNeeded(cgnsIn, cgnsOut)

contains
	!------------------------------------------------------------------------------------------
	!> �K�v�ȏꍇ�̂݁A�^�C���X�e�b�v�ʂ̊i�q���R�s�[����
	!> �^�C���X�e�b�v���ƂɊi�q���ω����Ȃ��v�Z���ʂ������ꍇ�͉������Ȃ�
	!------------------------------------------------------------------------------------------
	subroutine gridCopyIfNeeded(cgnsIn, cgnsOut)
		integer, intent(in) :: cgnsIn, cgnsOut
		integer :: zoneCopied, gridCopied
		integer :: baseId, zoneId, gridId;

		baseId = 1
		zoneId = 1
		gridId = 1

		! �܂��A�o�͐� (cgnsOut) �ɁAZone �����邩���ׂ�B�Ȃ���΍쐬����B
		call copyZoneIfNeeded(cgnsIn, cgnsOut, baseId, zoneId, zoneCopied)

		if (zoneCopied .eq. 0 ) then
			! Zone �̃R�s�[���s���Ȃ������B�i�q�͊��ɂ���B����ĉ������Ȃ��B
			return
		end if
		
		! �i�q���R�s�[���� (���W�l�̂�)
		call copyGrid(cgnsIn, cgnsOut, baseId, zoneId, gridId, gridId, gridCopied)
		! ��\���i�q�̍\�����R�s�[����
		call copyAllElements(cgnsIn, cgnsOut, baseId, zoneId)
		! �R�s�[���č쐬�����i�q�� iriclib �ŗ��p�ł���悤������
		call cg_iric_initgrid_mul_f(cgnsOut, zoneId, ier)

	end subroutine

	!------------------------------------------------------------------------------------------
	!> CI �̒l��ǂݍ��ށB�ǂݍ��߂��� ciIsSet �� 1 �ƂȂ�
	!------------------------------------------------------------------------------------------
	subroutine loadCI(cgnsOut, ciIsSet, ci)
		integer, intent(in) :: cgnsOut
		integer, intent(out) :: ciIsSet
		double precision, dimension(:), allocatable, intent(out) :: ci
		integer :: ier

		allocate(ci(gridSize))

		call cg_iric_read_grid_real_node_mul_f(cgnsOut, 'ChannelIndex', ci, ier)
		if (ier == 0) then
			ciIsSet = 1
		else
			ciIsSet = 0
		end if

	end subroutine

end subroutine

!------------------------------------------------------------------------------------------
!> �K�v�ȏꍇ�̂݁A�^�C���X�e�b�v�ʂ̊i�q���R�s�[����
!> �^�C���X�e�b�v���ƂɊi�q���ω����Ȃ��v�Z���ʂ������ꍇ�͉������Ȃ�
!------------------------------------------------------------------------------------------
subroutine gridCopyResultIfNeeded(solId, cgnsIn, cgnsOut)
	use grid

	implicit none

	integer, intent(in) :: solId
	integer, intent(in) :: cgnsIn, cgnsOut
	integer :: baseId, zoneId, gridId;
	integer :: gridCopied

	baseId = 1
	zoneId = 1
	gridId = solId + 1

	! �i�q���R�s�[����B1�Ԗڂ̊i�q�����͊i�q�Ȃ̂ŁA solId + 1 �Ԗڂ̊i�q���R�s�[����B
	call copyGrid(cgnsIn, cgnsOut, baseId, zoneId, gridId, gridId, gridCopied)

end subroutine

!------------------------------------------------------------------------------------------
!> �K�v�ȏꍇ�̂݁A�^�C���X�e�b�v�ʂ̊i�q���R�s�[����B (t = 999999 �p)
!> �^�C���X�e�b�v���ƂɊi�q���ω����Ȃ��v�Z���ʂ������ꍇ�͉������Ȃ�
!------------------------------------------------------------------------------------------
subroutine gridCopyResultLastIfNeeded(solCount, cgnsIn, cgnsOut)
	use common_params
	use grid

	implicit none

	integer, intent(in) :: solCount
	integer, intent(in) :: cgnsIn, cgnsOut
	integer :: baseId, zoneId, gridId;
	integer :: gridCopied
	character(len = strmax) newGridName

	baseId = 1
	zoneId = 1

	! �i�q���R�s�[����B1�Ԗڂ̊i�q�����͊i�q�Ȃ̂ŁA solId + 1 �Ԗڂ̊i�q���R�s�[����B
	call copyGridAs(cgnsIn, cgnsOut, baseId, zoneId, solCount + 1, solCount + 2, 'GridCoordinatesForLastSolution', gridCopied)

end subroutine

!------------------------------------------------------------------------------------------
!> �i�q�̖��O�̃��X�g����������
!------------------------------------------------------------------------------------------
subroutine gridWriteGridNames(cgnsOut)
	use common_params
	use grid
	
	integer, intent(in) :: cgnsOut
	integer :: baseId
	integer :: zoneId
	integer :: numgrids
	character(len = strmax), dimension(:), allocatable :: gridNames
	integer :: i, dim, ier
	integer, dimension(2) :: dims

	baseId = 1
	zoneId = 1

	call cg_ngrids_f(cgnsOut, baseId, zoneId, numGrids, ier)
	if (numGrids == 1) then
		return
	end if

	allocate(gridNames(numGrids - 1))
	do i = 1, numGrids - 1
		call cg_grid_read_f(cgnsOut, baseId, zoneId, i + 1, gridNames(i), ier)
	end do
	dim = 2
	dims(1) = strmax
	dims(2) = numgrids - 1
	call cg_goto_f(cgnsOut, baseId, ier, 'Zone_t', 1, 'ZoneIterativeData_t', 1, 'end')
	call cg_array_write_f('GridCoordinatesPointers', Character, dim, dims, gridNames, ier)

end subroutine

!------------------------------------------------------------------------------------------
!> �m�ۂ������������J��
!------------------------------------------------------------------------------------------
subroutine gridFree()
	use grid

	implicit none

	deallocate(roughness)
	deallocate(ci)
end subroutine
