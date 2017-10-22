program Euler
   implicit none

   double precision, allocatable :: x(:,:), xdot(:), Sol(:, :), Flux(:, :)
   double precision :: Time, OutputTime, DeltaT, DeltaX, CFL, ReportTime
   real :: StartTime, EndTime
   integer :: CurrentReportStep, NoOfReportSteps
   integer :: NumberOfNodes, InitialConditionChoice
   integer :: BoundaryConditionLeft, BoundaryConditionRight
   integer :: MovingMeshChoice

   ! Start timing procedure
   call cpu_time(StartTime)

   call ReadVariables(NumberOfNodes, NoOfReportSteps, CFL, InitialConditionChoice, MovingMeshChoice)

   allocate (x(NumberOfNodes, 2), xdot(NumberOfNodes))
   allocate (Sol(0:NumberOfNodes, 3), Flux(NumberOfNodes, 3))

   Time = 0.0d0
   CurrentReportStep = 1

   ! Create mesh
   call InitialMesh(x, xdot, NumberOfNodes, DeltaX)

   ! Set the initial condition
   call InitialCondition(x, Sol, NumberOfNodes, InitialConditionChoice, OutputTime &
                         , BoundaryConditionLeft, BoundaryConditionRight)

   ! Write time-stepping banner
   write (6, *)
   write (6, *) ' Advancing solution to time ', OutputTime
   write (6, *)
   write (6, *) '---------------------------------------------------------------------'
   write (6, *) '         Time                Time-step                               '
   write (6, *) '---------------------------------------------------------------------'

   open(unit=20,file='x_traj.m')
   write(20,*) x(:,1), Time   

   do while (Time < OutputTime)

      ! Calculate the time-step
      call CalculateTimeStep(DeltaT, CFL, x, Sol, NumberOfNodes)

      ! Adjust the time-step to reach report time
      ReportTime = dble(CurrentReportStep)*OutputTime/dble(NoOfReportSteps)
      DeltaT = min(Time + DeltaT, ReportTime) - Time

      ! Move the mesh
      if (MovingMeshChoice .gt. 1) then
         call LagrangianVelocity(Sol,x(:,1),xdot,NumberOfNodes)
      endif
      x(:,2) = x(:,1) + DeltaT * xdot

      ! Set the ghost cells
      call SetGhostCells(Sol, NumberOfNodes, BoundaryConditionLeft, BoundaryConditionRight)

      ! Calculate the flows
      call CalculateFlux(Sol, Flux, NumberOfNodes)

      ! Update the solution
      call UpdateSolution(x, Sol, Flux, NumberOfNodes, DeltaT)

      ! Increment the simulation time
      Time = Time + DeltaT

      ! Write out the current time and time-step size
      write (6, '(2f20.8)') Time, DeltaT

      ! Write the current solution
      if (abs(Time - ReportTime) < 1d-10) then
         call WriteSolution(x, Sol, NumberOfNodes)
         CurrentReportStep = CurrentReportStep + 1
         write(20,*) x(:,1), Time   
      endif
   enddo

   close(10)

   call WriteVariables(OutputTime, NoOfReportSteps)

   ! End timing procedure
   call cpu_time(EndTime)

   write (6, *)
   write (6, *) '---------------------------------------------------------------------'
   write (6, *) ' Computational time is ', EndTime - StartTime
   write (6, *) '---------------------------------------------------------------------'
   write (6, *)

end program Euler

subroutine ReadVariables(NumberOfNodes, NoOfReportSteps, CFL, InitialConditionChoice, MovingMeshChoice)
   implicit none
!-------------------------------------------------------------------------------
   double precision, intent(OUT) :: CFL
   integer, intent(OUT) :: NumberOfNodes, NoOfReportSteps, InitialConditionChoice
   integer, intent(OUT) :: MovingMeshChoice
!-------------------------------------------------------------------------------

   ! Read variables into program.
   open (unit=10, file='variables.data', status='old', form='formatted')
   read (10, *) NumberOfNodes
   read (10, *) NoOfReportSteps
   read (10, *) CFL
   read (10, *) InitialConditionChoice
   read (10, *) MovingMeshChoice
   close (10)

end subroutine ReadVariables

subroutine WriteVariables(OutputTime, NoOfReportSteps)

   implicit none
!------------------------------------------------------------------------------
   double precision, intent(IN) :: OutputTime
   integer, intent(IN) :: NoOfReportSteps
!------------------------------------------------------------------------------

   ! Write the solution variables to file
   open (unit=10, file='variables.m')
   write (10, *) OutputTime
   write (10, *) NoOfReportSteps
   close (10)

   return

end subroutine WriteVariables

subroutine InitialMesh(x, xdot, NumberOfNodes, DeltaX)
   implicit none
!-------------------------------------------------------------------------------
   integer, intent(in) :: NumberOfNodes
   double precision, intent(out) :: x(NumberOfNodes,2), xdot(NumberOfNodes), DeltaX
   double precision :: XMin, XMax
   integer :: i
!-------------------------------------------------------------------------------
   XMin = 0.0d0
   XMax = 1.0d0
   DeltaX = (XMax - XMin)/dble(NumberOfNodes - 1)

   do i = 1, NumberOfNodes
      x(i,1) = XMin + dble(i - 1)*DeltaX
      x(i,2) = x(i,1)
      xdot(i) = 0.0
   enddo

   return
end subroutine InitialMesh

subroutine InitialCondition(x, Sol, NumberOfNodes, InitialConditionChoice, OutputTime &
                            , BoundaryConditionLeft, BoundaryConditionRight)
   implicit none
!-------------------------------------------------------------------------------
   integer, intent(in) :: NumberOfNodes, InitialConditionChoice
   integer, intent(out) :: BoundaryConditionLeft, BoundaryConditionRight
   double precision, intent(out) :: OutputTime
   double precision, intent(out) :: x(NumberOfNodes,2)
   double precision, intent(out) :: Sol(0:NumberOfNodes, 3)
   integer :: i
!-------------------------------------------------------------------------------
   do i = 1, NumberOfNodes
      if (InitialConditionChoice == 1) then
         OutputTime = 0.1d0
         BoundaryConditionLeft = 1
         BoundaryConditionRight = 1
         ! Sod Problem
         if (x(i,1) < 0.5d0) then
            Sol(i, 1) = 1.0d0
            Sol(i, 2) = 0.0d0
            Sol(i, 3) = 2.5d0
         else
            Sol(i, 1) = 0.125d0
            Sol(i, 2) = 0.0d0
            Sol(i, 3) = 0.25d0
         endif
      elseif (InitialConditionChoice == 2) then
         ! Blast wave problem
         OutputTime = 0.038d0
         BoundaryConditionLeft = 2
         BoundaryConditionRight = 2
         if (x(i,1) <= 0.1d0) then
            Sol(i, 1) = 1.0d0
            Sol(i, 2) = 0.0d0
            Sol(i, 3) = 1000.0d0/0.4d0
         elseif (x(i,1) > 0.9) then
            Sol(i, 1) = 1.0d0
            Sol(i, 2) = 0.0d0
            Sol(i, 3) = 100.0d0/0.4d0
         else
            Sol(i, 1) = 1.0d0
            Sol(i, 2) = 0.0d0
            Sol(i, 3) = 0.01d0/0.4d0
         endif
      elseif (InitialConditionChoice == 3) then
         ! Lax problem
         OutputTime = 0.14d0
         BoundaryConditionLeft = 1
         BoundaryConditionRight = 1
         if (x(i,1) < 0.5d0) then
            Sol(i, 1) = 0.445d0
            Sol(i, 2) = 0.445d0*0.698d0
            Sol(i, 3) = 3.528d0/0.4d0 + 0.5d0*0.445d0*0.698d0*0.698d0
         else
            Sol(i, 1) = 0.5d0
            Sol(i, 2) = 0.0d0
            Sol(i, 3) = 0.571d0/0.4d0
         endif
      elseif (InitialConditionChoice == 4) then
         ! Modified Sod problem
         OutputTime = 0.2d0
         BoundaryConditionLeft = 1
         BoundaryConditionRight = 1
         if (x(i,1) < 0.5d0) then
            Sol(i, 1) = 1.0d0
            Sol(i, 2) = 1.0d0*0.75d0
            Sol(i, 3) = 1.0d0/0.4d0 + 0.5d0*1.0d0*0.75d0*0.75d0
         else
            Sol(i, 1) = 0.125d0
            Sol(i, 2) = 0.0d0
            Sol(i, 3) = 0.1d0/0.4d0
         endif
      elseif (InitialConditionChoice == 5) then
         ! 123 problem
         OutputTime = 0.15d0
         BoundaryConditionLeft = 1
         BoundaryConditionRight = 1
         if (x(i,1) < 0.5d0) then
            Sol(i, 1) = 1.0d0
            Sol(i, 2) = -2.0d0
            Sol(i, 3) = 0.4d0/0.4d0 + 0.5d0*1.0d0*2d0*2d0
         else
            Sol(i, 1) = 1.0d0
            Sol(i, 2) = 2.0d0
            Sol(i, 3) = 0.4d0/0.4d0 + 0.5d0*1.0d0*2d0*2d0
         endif
      endif
   enddo

   return
end subroutine InitialCondition

subroutine CalculateTimeStep(DeltaT, CFL, x, Sol, NumberOfNodes)
   implicit none
!-------------------------------------------------------------------------------
   integer, intent(in) :: NumberOfNodes
   double precision, intent(out) :: DeltaT
   double precision, intent(in) :: CFL, Sol(0:NumberOfNodes, 3), x(NumberOfNodes,2)
   double precision :: v, Lambda(3), a, Pressure, State(3)
   double precision, external :: Eos
   integer :: i, j
!-------------------------------------------------------------------------------

   v = 0.0d0
   do i = 1, NumberOfNodes - 1
      State = Sol(i, :)
      Pressure = Eos(State)
      a = sqrt(1.4d0*Pressure/Sol(i, 1))
      ! Get the Eigen-values
      Lambda(1) = (Sol(i, 2)/Sol(i, 1)) - a
      Lambda(2) = (Sol(i, 2)/Sol(i, 1))
      Lambda(3) = (Sol(i, 2)/Sol(i, 1)) + a
      do j = 1, 3
         v = max(v, abs(Lambda(j))/(x(i + 1,1) - x(i,1)))
      enddo
   enddo

!  v * (dt / dx) < CFL
   DeltaT = CFL/v

   return
end subroutine CalculateTimeStep

subroutine SetGhostCells(Sol, NumberOfNodes, BoundaryConditionLeft, BoundaryConditionRight)
   implicit none
!-------------------------------------------------------------------------------
   integer, intent(in) :: NumberOfNodes, BoundaryConditionLeft, BoundaryConditionRight
   double precision, intent(inout) :: Sol(0:NumberOfNodes, 3)
!-------------------------------------------------------------------------------
   ! Set left hand boundary
   if (BoundaryConditionLeft == 1) then
      Sol(0, :) = Sol(1, :)
   elseif (BoundaryConditionLeft == 2) then
      Sol(0, :) = Sol(1, :)
      Sol(0, 2) = -Sol(1, 2)
   endif

   ! Set right hand boundary
   if (BoundaryConditionRight == 1) then
      Sol(NumberOfNodes, :) = Sol(NumberOfNodes - 1, :)
   elseif (BoundaryConditionRight == 2) then
      Sol(NumberOfNodes, :) = Sol(NumberOfNodes - 1, :)
      Sol(NumberOfNodes, 2) = -Sol(NumberOfNodes - 1, 2)
   endif

   return
end subroutine SetGhostCells

subroutine CalculateFlux(Sol, Flux, NumberOfNodes)
   implicit none
!-------------------------------------------------------------------------------
   integer, intent(in) :: NumberOfNodes
   double precision, intent(out) :: Flux(NumberOfNodes, 3)
   double precision, intent(in) :: Sol(0:NumberOfNodes, 3)
   double precision :: SolLeft(3), SolRight(3), FluxL(3)
   integer :: i
!-------------------------------------------------------------------------------

   do i = 0, NumberOfNodes - 1
      SolLeft = Sol(i, :)
      SolRight = Sol(i + 1, :)
      call RoeFlux(FluxL, SolLeft, SolRight)
      Flux(i + 1, :) = FluxL
   enddo

   return
end subroutine CalculateFlux

subroutine FluxFunction(State, Flux)
   implicit none
!-------------------------------------------------------------------------------
   double precision, intent(in)  :: State(3)
   double precision, intent(out) :: Flux(3)
   double precision, external :: Eos
!-------------------------------------------------------------------------------
   Flux(1) = State(2)
   Flux(2) = State(2)*State(2)/State(1) + Eos(State)
   Flux(3) = (State(2)/State(1))*(State(3) + Eos(State))

end subroutine FluxFunction

function Eos(State) Result(Pressure)
   implicit none
!-------------------------------------------------------------------------------
   double precision, intent(in) :: State(*)
   double precision :: Pressure
   double precision :: Gamma = 1.4d0
!-------------------------------------------------------------------------------

! Relationships:
! --------------
! Pressure:         P = (Gamma-1)*(E-0.5*Rho*u^2)
! Enthalpy:         H = (E+P)/Rho = h+0.5*u^2
! Speed:            c = Gamma*P/Rho  = (Gamma-1)*(H-0.5*u^2)
! Internal energy:  e = E/Rho

   Pressure = (Gamma - 1d0)*(State(3) - 0.5d0*State(2)*State(2)/State(1))

   return
end function Eos

subroutine UpdateSolution(x, Sol, Flux, NumberOfNodes, DeltaT)
   implicit none
!-------------------------------------------------------------------------------
   integer, intent(in) :: NumberOfNodes
   double precision, intent(inout) :: x(NumberOfNodes,2)
   double precision, intent(in) :: Flux(NumberOfNodes, 3)
   double precision, intent(in) :: DeltaT
   double precision, intent(inout) :: Sol(0:NumberOfNodes, 3)
   double precision :: NewDeltaX, OldDeltaX
   integer :: i
!-------------------------------------------------------------------------------

   do i = 1, NumberOfNodes - 1
      OldDeltaX = x(i + 1,1) - x(i,1)
      NewDeltaX = x(i + 1,2) - x(i,2)
      Sol(i, :) = (OldDeltaX/NewDeltaX)*Sol(i, :) - (DeltaT/NewDeltaX)*(Flux(i + 1, :) - Flux(i, :))
   enddo

   x(:,1) = x(:,2)
   x(:,2) = 0.0d0;

   return
end subroutine UpdateSolution

subroutine WriteSolution(x, Sol, NumberOfNodes)
   implicit none
!-------------------------------------------------------------------------------
   integer, intent(in) :: NumberOfNodes
   double precision, intent(in) :: x(NumberOfNodes,2), Sol(0:NumberOfNodes, 3)
!-------------------------------------------------------------------------------
   ! Variables for writing output files
   character(LEN=40) :: Filename
   character(LEN=10) :: numbers = '0123456789'
   integer, parameter :: NoOfSF = 3
   integer :: TestNumber, i, units(NoOfSF)
   integer, save :: ReportStep = 1
   double precision, external :: Eos
!-------------------------------------------------------------------------------
   ! Construct the solution output filename
   TestNumber = ReportStep
   Filename = 'Solution'

   do i = 1, NoOfSF
      Units(i) = TestNumber/(10**(NoOfSF - i))
      TestNumber = TestNumber - (10**(NoOfSF - i))*units(i)
      Filename = trim(Filename)//numbers(units(i) + 1:units(i) + 1)
   enddo
   Filename = trim(Filename)//'.m'
   open (unit=10, file=Filename)
   do i = 1, NumberOfNodes - 1
      write (10, '(99f20.8)') 0.5d0*(x(i + 1,1) + x(i,1)), Sol(i, 1), Sol(i, 2)/Sol(i, 1), Eos(Sol(i,:)), Sol(i,3)/Sol(i,1) 
   enddo
   close (unit=10)

   ReportStep = ReportStep + 1

   return
end subroutine WriteSolution

subroutine RoeFlux(Flux, SolLeft, SolRight)
   implicit none
!-------------------------------------------------------------------------------
   double precision, intent(out) :: Flux(3)
   double precision, intent(in) :: SolLeft(3), SolRight(3)
   double precision :: Gamma = 1.4d0
   double precision :: K(3, 3)
   double precision :: Lambda(3), Alpha(3), DeltaSol(3), RoeAverage(3)
   double precision :: FluxLeft(3), FluxRight(3), a
   integer :: i
!-------------------------------------------------------------------------------

   ! Get the Roe average state
   call RoeStates(RoeAverage, SolLeft, SolRight)

   ! Get the Eigen-system
   call EigenSystem(Lambda, K, RoeAverage)

   ! Get the difference in the solution
   DeltaSol = SolRight - SolLeft

   a = sqrt((Gamma - 1d0)*(RoeAverage(3) - 0.5d0*RoeAverage(2)*RoeAverage(2)))

   ! Calculate the wave strengths. Solution of system K.Alpha = DeltaSol
   Alpha(2) = ((Gamma - 1.0d0)/(a*a))*(DeltaSol(1)*(RoeAverage(3) - (RoeAverage(2)*RoeAverage(2))) &
                                       + RoeAverage(2)*DeltaSol(2) - DeltaSol(3))
   Alpha(1) = (1.0d0/(2.0d0*a))*(DeltaSol(1)*(RoeAverage(2) + a) - DeltaSol(2) - a*Alpha(2))
   Alpha(3) = DeltaSol(1) - (Alpha(1) + Alpha(2))

   ! Calculate the flux function
   call FluxFunction(SolLeft, FluxLeft)
   call FluxFunction(SolRight, FluxRight)

   Flux = 0.5d0*(FluxLeft + FluxRight)
   do i = 1, 3
      Flux = Flux - 0.5d0*Alpha(i)*dabs(Lambda(i))*K(:, i)
   enddo

   return

end subroutine RoeFlux

subroutine RoeStates(RoeAverage, StateLeft, StateRight)
   implicit none
!-------------------------------------------------------------------------------
   double precision, intent(in) :: StateLeft(3), StateRight(3)
   double precision, intent(out) :: RoeAverage(3)
   double precision :: EntropyLeft, EntropyRight
   double precision :: RhoSqrtLeft, RhoSqrtRight, Denom
   double precision, external :: Eos
!-------------------------------------------------------------------------------

   EntropyLeft = (StateLeft(3) + Eos(StateLeft))/StateLeft(1)
   EntropyRight = (StateRight(3) + Eos(StateRight))/StateRight(1)

   RhoSqrtLeft = sqrt(StateLeft(1))
   RhoSqrtRight = sqrt(StateRight(1))
   Denom = RhoSqrtLeft + RhoSqrtRight

   ! Rho, Speed(u), Enthalpy(H)
   RoeAverage(1) = RhoSqrtLeft*RhoSqrtRight
   RoeAverage(2) = (RhoSqrtLeft*(StateLeft(2)/StateLeft(1)) + RhoSqrtRight*(StateLeft(2)/StateLeft(1)))/Denom
   RoeAverage(3) = (RhoSqrtLeft*EntropyLeft + RhoSqrtRight*EntropyRight)/Denom

   return

end subroutine RoeStates

subroutine EigenSystem(Values, Vectors, State)
   implicit none
!-------------------------------------------------------------------------------
   double precision, intent(in)  :: State(3)
   double precision, intent(out) :: Values(3), Vectors(3, 3)
   double precision :: a, gamma = 1.4d0
!-------------------------------------------------------------------------------

   a = sqrt((Gamma - 1d0)*(State(3) - 0.5d0*State(2)*State(2)))

   Values(1) = State(2) - a
   Values(2) = State(2)
   Values(3) = State(2) + a

   Vectors(1, 1) = 1.0d0
   Vectors(1, 2) = 1.0d0
   Vectors(1, 3) = 1.0d0
   Vectors(2, 1) = State(2) - a
   Vectors(2, 2) = State(2)
   Vectors(2, 3) = State(2) + a
   Vectors(3, 1) = State(3) - State(2)*a
   Vectors(3, 2) = 0.5d0*State(2)*State(2)
   Vectors(3, 3) = State(3) + State(2)*a

   return
end subroutine EigenSystem

subroutine LagrangianVelocity(u,x,xdot,NumberOfNodes)
   use linear_solvers

   implicit none
!---------------------------------------------------------------------------------      
!        GLOBAL VARIABLES
   integer, intent(IN) :: NumberOfNodes
   double precision, intent(INOUT), dimension(0:NumberOfNodes, 3) :: u
   double precision, intent(IN), dimension(NumberOfNodes) :: x
   double precision, intent(INOUT), dimension(NumberOfNodes) :: xdot
!---------------------------------------------------------------------------------      
!        LOCAL VARIABLES
   double precision, dimension(NumberOfNodes,NumberOfNodes) :: K
   double precision, dimension(NumberOfNodes) :: rhs, psi
   integer :: j    
!---------------------------------------------------------------------------------
   K=0.0d0; xdot=0.0d0; psi=0.0d0

   K(1,1) = -0.5d0*((u(2,1)+u(1,1))/(x(2)-x(1))) -0.5d0*((u(1,1)+u(2,1))/(x(1)+x(2)))
   K(2,1) =  0.5d0*((u(2,1)+u(1,1))/(x(2)-x(1)))
   rhs(1) = 0.5d0*( u(2,2) + u(1,2) )
   
   do j=2, NumberOfNodes-1
      K(j-1,j) = 0.5d0*((u(j,1)+u(j-1,1))/(x(j)-x(j-1)))
      K(j,j) = -0.5d0*((u(j,1)+u(j-1,1))/(x(j)-x(j-1))) - 0.5d0*((u(j+1,1)+u(j,1))/(x(j+1)-x(j)))
      K(j+1,j) = 0.5d0*((u(j+1,1)+u(j,1))/(x(j+1)-x(j)))
      rhs(j) = 0.5d0*( u(j+1,2) - u(j-1,2) )
   end do
        
   K(NumberOfNodes-1,NumberOfNodes) = 0.5d0*((u(NumberOfNodes,1)+u(NumberOfNodes-1,1))/(x(NumberOfNodes)-x(NumberOfNodes-1)))
   K(NumberOfNodes,NumberOfNodes) = -0.5d0*((u(NumberOfNodes,1)+u(NumberOfNodes-1,1))/(x(NumberOfNodes)-x(NumberOfNodes-1)))
   rhs(NumberOfNodes) = -0.5d0*( u(NumberOfNodes,2) + u(NumberOfNodes-1,2) )
                   
   call gaussian_elimination(K,psi,rhs,NumberOfNodes)

   call VelocityRecovery(x,xdot,psi,NumberOfNodes)        

   return

end subroutine LagrangianVelocity

subroutine VelocityRecovery(x,xdot,psi,NumberOfNodes)
   use linear_solvers

   implicit none
!---------------------------------------------------------------------------------
!        GLOBAL VARIABLES
   integer, intent(IN) :: NumberOfNodes
   double precision, intent(INOUT), dimension(NumberOfNodes) :: xdot
   double precision, intent(IN), dimension(NumberOfNodes) :: x
   double precision, intent(IN), dimension(NumberOfNodes) :: psi
!---------------------------------------------------------------------------------
!        LOCAL VARIABLES
   double precision, dimension(NumberOfNodes,NumberOfNodes) :: Mass
   double precision, dimension(NumberOfNodes) :: f
   integer :: j
!--------------------------------------------------------------------------------
   Mass=0.0d0; f=0.0d0; xdot=0.0d0
                
   Mass(1,1) = (1.0d0/3.0d0)*(x(2)-x(1))
   Mass(2,1) = (1.0d0/6.0d0)*(x(2)-x(1))
   f(1) = 0.5d0*(psi(2)-psi(1))
   
   do j=2, NumberOfNodes-1
      Mass(j-1,j) = (1.0d0/6.0d0)*(x(j)-x(j-1))
      Mass(j,j) = (1.0d0/3.0d0)*(x(j+1)-x(j-1))
      Mass(j+1,j) = (1.0d0/6.0d0)*(x(j+1)-x(j))
      f(j) = 0.5d0*(psi(j+1)-psi(j-1))
   end do
        
   Mass(NumberOfNodes-1,NumberOfNodes) = (1.0d0/6.0d0)*(x(NumberOfNodes)-x(NumberOfNodes-1))
   Mass(NumberOfNodes,NumberOfNodes) = (1.0d0/3.0d0)*(x(NumberOfNodes)-x(NumberOfNodes-1))
   f(NumberOfNodes) = 0.5d0*(psi(NumberOfNodes)-psi(NumberOfNodes-1))
        
   call gaussian_elimination(Mass,xdot,f,NumberOfNodes)
        
   return

end subroutine VelocityRecovery
