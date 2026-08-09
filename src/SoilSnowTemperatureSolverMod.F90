module SoilSnowTemperatureSolverMod

!!! Compute soil and snow layer temperature using tri-diagonal matrix solution
!!! Dependent on the output from SoilSnowThermalDiffusion subroutine

  use Machine
  use NoahmpVarType
  use ConstantDefineMod
  use MatrixSolverTriDiagonalMod, only : MatrixSolverTriDiagonal

  implicit none

contains

  subroutine SoilSnowTemperatureSolver(noahmp, TimeStep, MatLeft1, MatLeft2, MatLeft3, MatRight)
#ifdef NOAHMP_ACC_COLUMNS
!$acc routine seq
#endif

! ------------------------ Code history --------------------------------------------------
! Original Noah-MP subroutine: HSTEP
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! ----------------------------------------------------------------------------------------

    implicit none

! in & out variables
    type(noahmp_type)     , intent(inout) :: noahmp
    real(kind=kind_noahmp), intent(in)    :: TimeStep                             ! timestep (may not be the same as model timestep)
#ifdef NOAHMP_ACC_COLUMNS
    ! Bounds must match the actual arrays declared in SoilSnowTemperatureMainMod.
    ! A plain dimension(:) dummy rebases to 1 and sends the negative snow-layer
    ! indices (LoopInd from NumSnowLayerNeg+1) off the front of the array.
    real(kind=kind_noahmp), dimension(-NoahmpAccMaxSnowLayer+1:NoahmpAccMaxSoilLayer), intent(inout) :: MatRight
    real(kind=kind_noahmp), dimension(-NoahmpAccMaxSnowLayer+1:NoahmpAccMaxSoilLayer), intent(inout) :: MatLeft1
    real(kind=kind_noahmp), dimension(-NoahmpAccMaxSnowLayer+1:NoahmpAccMaxSoilLayer), intent(inout) :: MatLeft2
    real(kind=kind_noahmp), dimension(-NoahmpAccMaxSnowLayer+1:NoahmpAccMaxSoilLayer), intent(inout) :: MatLeft3
#else
    real(kind=kind_noahmp), allocatable, dimension(:), intent(inout) :: MatRight  ! right-hand side term of the matrix
    real(kind=kind_noahmp), allocatable, dimension(:), intent(inout) :: MatLeft1  ! left-hand side term of the matrix
    real(kind=kind_noahmp), allocatable, dimension(:), intent(inout) :: MatLeft2  ! left-hand side term of the matrix
    real(kind=kind_noahmp), allocatable, dimension(:), intent(inout) :: MatLeft3  ! left-hand side term of the matrix
#endif

! local variable
    integer                                           :: LoopInd                  ! layer loop index 
#ifdef NOAHMP_ACC_COLUMNS
    real(kind=kind_noahmp), dimension(-NoahmpAccMaxSnowLayer+1:NoahmpAccMaxSoilLayer) :: MatRightTmp
    real(kind=kind_noahmp), dimension(-NoahmpAccMaxSnowLayer+1:NoahmpAccMaxSoilLayer) :: MatLeft3Tmp
#else
    real(kind=kind_noahmp), allocatable, dimension(:) :: MatRightTmp              ! temporary MatRight matrix coefficient
    real(kind=kind_noahmp), allocatable, dimension(:) :: MatLeft3Tmp              ! temporary MatLeft3 matrix coefficient
#endif

! --------------------------------------------------------------------
    associate(                                                                &
              NumSoilLayer        => noahmp%config%domain%NumSoilLayer       ,& ! in,    number of soil layers
              NumSnowLayerMax     => noahmp%config%domain%NumSnowLayerMax    ,& ! in,    maximum number of snow layers
              NumSnowLayerNeg     => noahmp%config%domain%NumSnowLayerNeg    ,& ! in,    actual number of snow layers (negative)
              TemperatureSoilSnow => noahmp%energy%state%TemperatureSoilSnow  & ! inout, snow and soil layer temperature [K]
             )
! ----------------------------------------------------------------------

    ! initialization
#ifndef NOAHMP_ACC_COLUMNS
    if (.not. allocated(MatRightTmp)) allocate(MatRightTmp(-NumSnowLayerMax+1:NumSoilLayer))
    if (.not. allocated(MatLeft3Tmp)) allocate(MatLeft3Tmp(-NumSnowLayerMax+1:NumSoilLayer))
#endif
    MatRightTmp = 0.0
    MatLeft3Tmp = 0.0

    ! update tri-diagonal matrix elements
    do LoopInd = NumSnowLayerNeg+1, NumSoilLayer
       MatRight(LoopInd) =       MatRight(LoopInd) * TimeStep
       MatLeft1(LoopInd) =       MatLeft1(LoopInd) * TimeStep
       MatLeft2(LoopInd) = 1.0 + MatLeft2(LoopInd) * TimeStep
       MatLeft3(LoopInd) =       MatLeft3(LoopInd) * TimeStep
    enddo

    ! copy values for input variables before call to rosr12
    do LoopInd = NumSnowLayerNeg+1, NumSoilLayer
       MatRightTmp(LoopInd) = MatRight(LoopInd)
       MatLeft3Tmp(LoopInd) = MatLeft3(LoopInd)
    enddo

    ! solve the tri-diagonal matrix equation
    call MatrixSolverTriDiagonal(MatLeft3,MatLeft1,MatLeft2,MatLeft3Tmp,MatRightTmp,&
                                 MatRight,NumSnowLayerNeg+1,NumSoilLayer,NumSnowLayerMax,&
                                 lbound(MatLeft3,1))

    ! update snow & soil temperature
    do LoopInd = NumSnowLayerNeg+1, NumSoilLayer
       TemperatureSoilSnow(LoopInd) = TemperatureSoilSnow(LoopInd) + MatLeft3(LoopInd)
    enddo

    ! deallocate local arrays to avoid memory leaks
#ifndef NOAHMP_ACC_COLUMNS
    deallocate(MatRightTmp)
    deallocate(MatLeft3Tmp)
#endif

    end associate

  end subroutine SoilSnowTemperatureSolver

end module SoilSnowTemperatureSolverMod
