!====================================================================================================
!
!> 格子の読み込みやコピーを行う処理を行うモジュール
!
!====================================================================================================

module grid
	use common_params

	implicit none

	!include 'cgnslib_f.h'
	!include 'iriclib_f.h'

	!> 格子点の数
	integer :: gridSize
	!> 粗度係数
	double precision, dimension(:), allocatable :: roughness
	!> ci が読み込まれていたら 1, いなかったら 0
	integer :: ciIsSet
	!> Channel Index (HSI解析用)
	double precision, dimension(:), allocatable :: ci

contains

	!------------------------------------------------------------------------------------------
	!> 格子の配列をコピー
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

		! 配列のメタ情報を取得
		call cg_goto_f(cgnsIn, baseId, ier, 'Zone_t', zoneId, 'GridCoordinates_t', gridId, 'end')
		call cg_array_info_f(arrayId, arrayName, dataType, numDims, dims, ier)

		dataSize = 1
		do i = 1, numDims
			dataSize = dataSize * dims(i)
		end do

		! 配列の値を読み込み
		allocate(arrayData(dataSize))
		call cg_array_read_f(arrayId, arrayData, ier)

		! 配列の値を書き込み
		call cg_goto_f(cgnsOut, baseId, ier, 'Zone_t', zoneId, 'GridCoordinates_t', gridIdTo, 'end')
		call cg_array_write_f(arrayName, dataType, numDims, dims, arrayData, ier)
		deallocate(arrayData)

	end subroutine

	!------------------------------------------------------------------------------------------
	!> 格子の全配列をコピー
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
	!> Zoneの Elements_t 要素をコピー
	!------------------------------------------------------------------------------------------
	subroutine copyElement(cgnsIn, cgnsOut, baseId, zoneId, elemId)
		integer, intent(in) :: cgnsIn, cgnsOut
		integer, intent(in) :: baseId, zoneId, elemId

		character(len = strmax) :: sectionName
		integer :: sectionType, start, end, nbndry, parent_flag, dataSize, parentData
		integer, dimension(:), allocatable :: data
		integer :: tmpElemId

		integer :: ier

		! section の情報を読み込み
		call cg_section_read_f(cgnsIn, baseId, zoneId, elemId, sectionName, &
			sectionType, start, end, nbndry, parent_flag, ier)
		call cg_elementdatasize_f(cgnsIn, baseId, zoneId, elemId, dataSize, ier)
		allocate(data(dataSize))
		call cg_elements_read_f(cgnsIn, baseId, zoneId, elemId, data, parentData, ier)

		! section の情報を書き込み
		call cg_section_write_f(cgnsOut, baseId, zoneId, sectionName, sectionType, &
			start, end, nbndry, data, tmpElemId, ier)

	end subroutine

	!------------------------------------------------------------------------------------------
	!> Zoneの全 Elements_t 要素をコピー
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
	!> 該当する Zone がなければ、コピーする。
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
			! Zone をコピーするべき時ではない
			zoneCopied = 0
			return
		end if

		! まず、コピー元の Zone の情報を読み込む
		call cg_zone_read_f(cgnsIn, baseId, zoneId, zoneName, zoneSize, ier)
		call cg_zone_type_f(cgnsIn, baseId, zoneId, zoneType, ier)

		! 同じ情報をコピー先に書き込む
		call cg_zone_write_f(cgnsOut, baseId, zoneName, zoneSize, zoneType, tmpZoneId, ier)
		zoneCopied = 1

	end subroutine

	!------------------------------------------------------------------------------------------
	!> 格子をコピーする。
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
			! コピーする格子がない
			gridCopied = 0
			return
		end if

		call cg_ngrids_f(cgnsOut, baseId, zoneId, numGrids, ier)
		
		if (numGrids /= gridIdTo - 1) then
			! 格子をコピーするべき時ではない
			gridCopied = 0
			return
		end if

		! GridCoordinates_t の情報を読み込む
		call cg_grid_read_f(cgnsIn, baseId, zoneId, gridId, coordName, ier)
		! GridCoordinates_t を作成
		call cg_grid_write_f(cgnsOut, baseId, zoneId, coordName, tmpGridId, ier)

		! GridCoordinate_t 以下の全ての 配列をコピー
		call copyAllArraysUnderGC(cgnsIn, cgnsOut, baseId, zoneId, gridId, gridIdTo)

		gridCopied = 1

	end subroutine

	!------------------------------------------------------------------------------------------
	!> 名前をつけて格子をコピーする。
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
			! コピーする格子がない
			gridCopied = 0
			return
		end if

		call cg_ngrids_f(cgnsOut, baseId, zoneId, numGrids, ier)
		
		if (numGrids /= gridIdTo - 1) then
			! 格子をコピーするべき時ではない
			gridCopied = 0
			return
		end if

		! GridCoordinates_t の情報を読み込む
		call cg_grid_read_f(cgnsIn, baseId, zoneId, gridId, coordName, ier)
		! GridCoordinates_t を作成
		call cg_grid_write_f(cgnsOut, baseId, zoneId, newName, tmpGridId, ier)

		! GridCoordinate_t 以下の全ての 配列をコピー
		call copyAllArraysUnderGC(cgnsIn, cgnsOut, baseId, zoneId, gridId, gridIdTo)

		gridCopied = 1

	end subroutine

end module

!------------------------------------------------------------------------------------------
!> 格子を読み込み、格子点の数を gridSize に出力する
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
	double precision, dimension(:,:), allocatable :: roughness_cell !< 粗度係数 (セル)
	double precision :: tmp_roughness

	! BASE, ZONE は基本的に1

	baseId = 1
	zoneId = 1

	! 格子情報を読み込み
	call cg_zone_read_f(cgnsIn, baseId, zoneId, zoneName, zoneSize, ier)
	call cg_zone_type_f(cgnsIn, baseId, zoneId, zoneType, ier)

	! 格子点の数を gridSize に格納
	if (zoneType .eq. 2) then
		gridSize = zoneSize(1) * zoneSize(2)
	else 
		gridSize = zoneSize(1)
	end if

	! Roughness を読み込む
	allocate(roughness(gridSize))
	roughness = 0

	if (zoneType .eq. 2) then
		! まず、roughness_cell の読み込みを試みる。
		allocate(roughness_cell(zoneSize(1) - 1, zoneSize(2) - 1))
		call cg_iric_read_grid_real_node_mul_f(cgnsIn, 'roughness_cell', roughness_cell, ier)
		if (ier .eq. 0) then
			! 読み込みが成功した。セルで定義されているので、その格子点を囲む4つのセルでの
			! 粗度係数の最大値を採用する
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
			! 読み込みが失敗した。格子点で定義されていると仮定
			call cg_iric_read_grid_real_node_mul_f(cgnsIn, 'Roughness', roughness, ier)
			if (ier .ne. 0) then
				! 読み込みに失敗したら、 "roughness" も試す
				call cg_iric_read_grid_real_node_mul_f(cgnsIn, 'roughness', roughness, ier)
			end if
		end if
	else 
		! roughness は、格子点で定義されていると仮定
		call cg_iric_read_grid_real_node_mul_f(cgnsIn, 'Roughness', roughness, ier)
		if (ier .ne. 0) then
			! 読み込みに失敗したら、 "roughness" も試す
			call cg_iric_read_grid_real_node_mul_f(cgnsIn, 'roughness', roughness, ier)
		end if
		if (ier .ne. 0) then
			! 読み込みに失敗したら、 "initial_cd" も試す (SToRM 用)
			call cg_iric_read_grid_real_node_mul_f(cgnsIn, 'initial_cd', roughness, ier)
		end if
	end if

	! CI の値を読み込む
	call loadCI(cgnsOut, ciIsSet, ci)

	! 最後に、必要があれば格子を cgnsIn から cgnsOut にコピーする
	call gridCopyIfNeeded(cgnsIn, cgnsOut)

contains
	!------------------------------------------------------------------------------------------
	!> 必要な場合のみ、タイムステップ別の格子をコピーする
	!> タイムステップごとに格子が変化しない計算結果だった場合は何もしない
	!------------------------------------------------------------------------------------------
	subroutine gridCopyIfNeeded(cgnsIn, cgnsOut)
		integer, intent(in) :: cgnsIn, cgnsOut
		integer :: zoneCopied, gridCopied
		integer :: baseId, zoneId, gridId;

		baseId = 1
		zoneId = 1
		gridId = 1

		! まず、出力先 (cgnsOut) に、Zone があるか調べる。なければ作成する。
		call copyZoneIfNeeded(cgnsIn, cgnsOut, baseId, zoneId, zoneCopied)

		if (zoneCopied .eq. 0 ) then
			! Zone のコピーが行われなかった。格子は既にある。よって何もしない。
			return
		end if
		
		! 格子をコピーする (座標値のみ)
		call copyGrid(cgnsIn, cgnsOut, baseId, zoneId, gridId, gridId, gridCopied)
		! 非構造格子の構造をコピーする
		call copyAllElements(cgnsIn, cgnsOut, baseId, zoneId)
		! コピーして作成した格子を iriclib で利用できるよう初期化
		call cg_iric_initgrid_mul_f(cgnsOut, zoneId, ier)

	end subroutine

	!------------------------------------------------------------------------------------------
	!> CI の値を読み込む。読み込めたら ciIsSet は 1 となる
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
!> 必要な場合のみ、タイムステップ別の格子をコピーする
!> タイムステップごとに格子が変化しない計算結果だった場合は何もしない
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

	! 格子をコピーする。1番目の格子が入力格子なので、 solId + 1 番目の格子をコピーする。
	call copyGrid(cgnsIn, cgnsOut, baseId, zoneId, gridId, gridId, gridCopied)

end subroutine

!------------------------------------------------------------------------------------------
!> 必要な場合のみ、タイムステップ別の格子をコピーする。 (t = 999999 用)
!> タイムステップごとに格子が変化しない計算結果だった場合は何もしない
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

	! 格子をコピーする。1番目の格子が入力格子なので、 solId + 1 番目の格子をコピーする。
	call copyGridAs(cgnsIn, cgnsOut, baseId, zoneId, solCount + 1, solCount + 2, 'GridCoordinatesForLastSolution', gridCopied)

end subroutine

!------------------------------------------------------------------------------------------
!> 格子の名前のリストを書き込む
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
!> 確保したメモリを開放
!------------------------------------------------------------------------------------------
subroutine gridFree()
	use grid

	implicit none

	deallocate(roughness)
	deallocate(ci)
end subroutine
