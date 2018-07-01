      PROGRAM RANDOMASSIGNMENT

!------------------DESCRIPTION-------------------------------------
!     This program performs two problems involving random numbers:
!     1) It calculates sets of k random numbers, where k = 2, 4, 8, 16, and 32.
!     For k = 32, the mean and std are calculated, and the probability density is
!     compared to the normal distribution.
!
!     Inputs: NONE
!     Outputs: See above Descriptions.
!
!     2) Determines the value of pi based on a drawing of random points
!     from a box surrounding a fourth of a circle. The ratio of points inside
!     the fourth circle versus in the box overall should roughly be equal to
!     the ratio of the areas of the circle and box.

!     Inputs: NONE
!     Outputs: The value of Pi.

!------------------------------------------------------------------

      IMPLICIT NONE
      REAL, DIMENSION (2) :: randarray2    !2 random numbers
      REAL, DIMENSION (4) :: randarray4    !4 random numbers
      REAL, DIMENSION (8) :: randarray8    !8 random numbers
      REAL, DIMENSION (16) :: randarray16    !16 random numbers
      REAL, DIMENSION (32) :: randarray32    !32 random numbers
      INTEGER, DIMENSION (32) :: seed    !values for seeding the random numbers
      REAL :: randavg, randstd, normaldist    !the average and stddev of 32 random numbers
      INTEGER :: values(1:8), k    !used to make note of the date and time, for seeding
!     The following are used in part 2 of the assignment

      REAL, DIMENSION (100) :: randxarray102, randyarray102
      REAL, DIMENSION (10000) :: randxarray104, randyarray104
      REAL, DIMENSION (1000000) :: randxarray106, randyarray106
      REAL, DIMENSION (100000000) :: randxarray108, randyarray108
      REAL :: sumsqrd, pi, error, pierror
      INTEGER :: i, numincircle


      open(72, file='RandomAssignment.out', form = 'formatted')

      call date_and_time(values=values)    !get the date and time to use as the seed

      k = 32   ! highest number of random values

      call random_seed(size=k)    !set up the random seed from the date and time values
      seed(1:32) = values(8)
      call random_seed(put=seed)

!----------PART ONE---------COMPARISON TO NORMAL DISTRIBUTION-----------------------


!     Generate the random number arrays
      call random_number(randarray2)
      call random_number(randarray4)
      call random_number(randarray8)
      call random_number(randarray16)
      call random_number(randarray32)

      write(72, *) 'Two Random Numbers:'
      write(72, '(F8.2,X)')  randarray2
      write(72, *) ''
      write(72, *) 'Four Random Numbers:'
      write(72, '(F8.2,X)')  randarray4
      write(72, *) ''
      write(72, *) 'Eight Random Numbers:' 
      write(72, '(F8.2,X)')  randarray8
      write(72, *) ''
      write(72, *) 'Sixteen Random Numbers:'
      write(72, '(F8.2,X)')  randarray16
      write(72, *) ''
      write(72, *) 'Thirty-Two Random Numbers:' 
      write(72, '(F8.2,X)')  randarray32
      write(72, *) ''

!     Find the average and standard deviation of the last set.
      randavg = sum(randarray32)/32.d0
      randstd = sqrt(sum((randarray32-randavg)**2)/32)

      write(72, *) 'The average of the thirty-two: ', randavg
      write(72, *) 'Standard Deviation: ', randstd


!-------------PART TWO------PI FINDER----------------------------------
      numincircle = 0   !Number of points inside the quarter circle
      sumsqrd = 0.d0    !The sum squared of the x and y coordinates
      pi = 0.d0    !The value of pi

      write(72, *) ''
      write(72, *) 'CALCULATE THE VALUE OF PI'

      call random_number (randxarray102)
      call random_number (randyarray102)

      do i = 1, 100
         ! x^2 + y^2 <= 1 to be inside the circle
         sumsqrd = (randxarray102(i)**2 + randyarray102(i)**2)
         if (sumsqrd.LE.1d0) numincircle = numincircle + 1 !Note if inside circle
      end do   

      pi = (numincircle/100.d0) * 4     !from pi*(1^2)/4 = number_in/total_number
      error = pierror(pi)

      write(72, *) 'NUM_POINTS   CALCULATED   ERROR'
      write(72, '(I12, 2(3X,f12.8))') 100, pi, error

      numincircle = 0  !reinitialize numincircle
      call random_number (randxarray104)
      call random_number (randyarray104)

      do i = 1, 10000
         ! x^2 + y^2 <= 1 to be inside the circle
         sumsqrd = (randxarray104(i)**2 + randyarray104(i)**2)
         if (sumsqrd.LE.1d0) numincircle = numincircle + 1 !Note if inside circle
      end do   

      pi = (numincircle/10000.d0) * 4     !from pi*(1^2)/4 = number_in/total_number
      error = pierror(pi)

      write(72, '(I12, 2(3X,f12.8))') 10000, pi, error

      numincircle = 0  !reinitialize numincircle
      call random_number (randxarray106)
      call random_number (randyarray106)

      do i = 1, 1000000
         ! x^2 + y^2 <= 1 to be inside the circle
         sumsqrd = (randxarray106(i)**2 + randyarray106(i)**2)
         if (sumsqrd.LE.1d0) numincircle = numincircle + 1 !Note if inside circle
      end do   

      pi = (numincircle/1000000.d0) * 4     !from pi*(1^2)/4 = number_in/total_number
      error = pierror(pi)

      write(72, '(I12, 2(3X,f12.8))') 1000000, pi, error

      numincircle = 0  !reinitialize numincircle
      call random_number (randxarray108)
      call random_number (randyarray108)

      do i = 1, 100000000
         ! x^2 + y^2 <= 1 to be inside the circle
         sumsqrd = (randxarray108(i)**2 + randyarray108(i)**2)
         if (sumsqrd.LE.1d0) numincircle = numincircle + 1 !Note if inside circle
      end do   

      pi = (numincircle/100000000.d0) * 4     !from pi*(1^2)/4 = number_in/total_number
      error = pierror(pi)

      write(72, '(I12, 2(3X,f12.8))') 100000000, pi, error


      END PROGRAM

      function normaldist(x, avg, stdev)
      IMPLICIT NONE
      REAL,  INTENT(IN) :: x, avg, stdev
      REAL :: normaldist

      normaldist = (1d0/sqrt(6.28*stdev**2))*
     2 exp(-((x-avg)**2)/(2d0*stdev**2))
      end function

      function pierror(calculated)
      IMPLICIT NONE
      REAL, INTENT(IN) :: calculated   !The calculated value of pi
      REAL :: pierror   !The percent different between the calculated and actual
      REAL, PARAMETER :: actual = 4*ATAN(1.d0)   !The actual value of pi
      pierror = 100*ABS(calculated - actual)/actual
      end function
