      PROGRAM ODE

      IMPLICIT NONE
      REAL, DIMENSION (11) :: x, y, y_an
      REAL, DIMENSION (21) :: x2, y2
      REAL, DIMENSION (51) :: x3, y3
      REAL, DIMENSION (101) :: x4, y4
      REAL, DIMENSION (201) :: x5, y5
      INTEGER :: i
      REAL :: fxn, analyticalsln, difference, diff 

      open(72, file='ODE_Assignment1.out', form='formatted')
      write(72,*) 'x_n   y_n   y_n_Analytical   Difference'

      x(1:11) = (/ 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,
     2 0.8, 0.9, 1.0 /)     !x is the input
      y(1) = 1d0    !y is the calculated solution
      y_an(1) = 1d0   !The analytic solution
      diff = 0      !The percent difference between actual and calculated

      write(72, '(4F8.2,X)') x(1), y(1), y_an(1), diff

!Step size 0.1
      do i = 1,10
         y(i+1) = fxn(x(i),y(i))*0.1 + y(i)   !Calculate next step
         y_an(i+1) = analyticalsln(x(i+1))    !Anal of next step
         diff = 100*(abs(y(i+1)-y_an(i+1))/(y_an(i+1)))
         write(72, '(4F8.2,X)') x(i+1), y(i+1), y_an(i+1), diff
      end do 

 
      write(72,*) ''
      write(72,*) 'Comparison of Different Stepsizes:'
      write(72,*) 'StepSize  0.1_Diff  1.0_Diff'


!Step size 0.05
      x2(1) = 0d0
      y2(1) = 1d0

      do i = 1,21
         x2(i+1) = i*0.05
         y2(i+1) = fxn(x2(i),y2(i))*0.05 + y2(i) 
      end do
      write(72, '(F3.2, 4F8.2)') 0.05, difference(y2(3),y_an(2)), 
     2 difference(y2(21),y_an(11))

!Step size 0.02
      x3(1) = 0d0
      y3(1) = 1d0

      do i = 1,51
         x3(i+1) = i*0.02
         y3(i+1) = fxn(x3(i),y3(i))*0.02 + y3(i)
      end do
      write(72, '(F3.2, 4F8.2)') 0.02, difference(y3(6),y_an(2)), 
     2 difference(y3(51),y_an(11))

!Step size 0.01
      x4(1) = 0d0
      y4(1) = 1d0

      do i = 1,101
         x4(i+1) = i*0.01
         y4(i+1) = fxn(x4(i),y4(i))*0.01 + y4(i)
      end do
      write(72, '(F3.2, 4F8.2)') 0.01, difference(y4(11),y_an(2)), 
     2 difference(y4(101),y_an(11))

!Step size 0.005
      x5(1) = 0d0
      y5(1) = 1d0

      do i = 1,201
         x5(i+1) = i*0.005
         y5(i+1) = fxn(x5(i),y5(i))*0.005 + y5(i)
      end do
      write(72, '(F4.3, 4F8.2)') 0.005, difference(y5(21),y_an(2)), 
     2 difference(y5(201),y_an(11))

      END PROGRAM

      function fxn(x, y)
      IMPLICIT NONE
      REAL, INTENT(IN) :: x, y
      REAL :: fxn
!The function we are solving
      fxn = 2d0*(y**2 + 1d0)/(x**2 + 4d0)
      end function


      function analyticalsln(x)
      IMPLICIT NONE
      REAL, INTENT(IN) :: x
      REAL :: analyticalsln
!Calculates the analytic solution
      analyticalsln = (2 + x)/(2 - x)
      end function


      function difference(y, y_an)
      IMPLICIT NONE
      REAL, INTENT(IN) :: y, y_an
      REAL :: difference
!Calculates the oercent difference between y and y_an
      difference = 100*(abs(y-y_an)/y_an)
      end function
