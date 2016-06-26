program dusty   
implicit none 
 
    real (kind=8) :: seed

    interface
        real (kind=8) function conrand(seed)
          real (kind=8) :: seed
        end function conrand    
    end interface  
  
      integer MAXDIM
      parameter (MAXDIM = 50)   

      integer IA(MAXDIM), N, i, j, k, ival
      double precision AV(MAXDIM), BV(MAXDIM), CV(MAXDIM)
      double precision OP(MAXDIM, MAXDIM), ID(MAXDIM, MAXDIM)
      double precision AM(MAXDIM, MAXDIM), BM(MAXDIM, MAXDIM)
      double precision CM(MAXDIM, MAXDIM), DM(MAXDIM, MAXDIM)
      double precision check, BOT, TOP, HOLDA, HOLDB, TRACE3
      double precision sub, jn, jn1, jn2
	  
#ifdef NOPRECISIONLOSS
      double precision start, finish, sum
#else
      real start, finish, sum
#endif

      double precision trig
      external trig

      N = MAXDIM
      seed = 1.0
      call cpu_time(start)
      !call srand(1)

!     Fill arrays

! Loop 10 Series -- Filling Arrays

      do i = 1, N
        jn1 = conrand(seed)
        jn2 = conrand(seed)
        jn = jn1 * (-1) ** (mod(int(10 * jn2), N))
		CV(i) = 0
        AV(i) = bessel_jn(0, jn)
      end do 
	  
      do i = 1, N
        jn1 = conrand(seed)
        jn2 = conrand(seed)
        jn = jn1 * (-1) ** (mod(int(10 * jn2), N))
        BV(i) = bessel_jn(1, jn)
      end do
     
      check = 0.0

      do i = 1, N
        ival = N
        check = check + AV(i) * BV(i)
        call idcheck(ival, check, AV, BV, ID)
      end do

! Compute |AV><BV|
#ifdef STRIDEACCESS
    do j = 1, N
        do i = 1, N
          call idcheck(N,check,AV,BV,ID)
          if (check > 0.5) then
             OP(i,j) = AV(i) * BV(j) / BV(i)
          else
             OP(i,j) = AV(j) * BV(i) / BV(j) 
          endif
       end do
	   !print *, i
       IA(i) = i
    end do
#else
    do i = 1, N
        do j = 1, N
          call idcheck(N,check,AV,BV,ID)
          if (check > 0.5) then
             OP(i,j) = AV(i) * BV(j) / BV(i)
          else
             OP(i,j) = AV(j) * BV(i) / BV(j) 
          endif
       end do
	   !print *, i
       IA(i) = i
    end do
#endif

	!loop 10
		
    do i = 1, N
        do j = 0, i, 8
            IA(I) = mod(mod(i+j,N),N)+1 !eliminar um mod
			!print *, IA(I) 
        end do
    end do

! Loop 20
    do i = 1, N
        call idcheck(N,check,AV,BV,ID)
        CV(IA(I)) = (AV(IA(I)) + BV(IA(I))) / check
		!print *, CV(IA(I))
    end do


! Loop 30

    do i = 2, N
        call idcheck(N,check,AV,BV,ID)
        AV(i) = AV(i-1) * BV(i) + CV(i)
		!print *, AV(i)
    end do

! Loop 40
    do i = 1, N
        call idcheck(N,check,AV,BV,ID)
        do j = 1, N
           if (check > 0.5) then
              BOT = OP(i,j)
              TOP = AV(j) * BV(j)
              HOLDA = AV(j)
              AV(j) = BV(j) + CV(j) / (TOP-BOT) * ID(i,i)
              BV(j) = HOLDA + CV(j) / (TOP-BOT) * ID(j,j)
              AM(i,j) = AV(j) * trig(IA(i), IA(j))
              BM(i,j) = BV(j) * trig(IA(j), IA(i))
           else
              BOT = OP(i,j)
              TOP = AV(j) * BV(j)
              HOLDA = AV(j)
              AV(j) = BV(j) - CV(j) / (TOP-BOT) * ID(j,j)
              BV(j) = HOLDA - CV(j) / (TOP-BOT) * ID(i,i)
              AM(i,j) = AV(j) / trig(IA(i), IA(j))
              BM(i,j) = BV(j) / trig(IA(j), IA(i))
           endif
        end do
    end do

#ifdef FOUR
! Loop 50
      do i = 1, N
         do  j = i + 1, N
            CM(i, j) = 0.0
            do k = 1, N
                CM(i, j) = CM(i, j) - AM(i, k) * BM(k, j) / check
            end do
         end do
      end do


      do i = 1, N
         do  j = 1, i
            CM(i, j) = 0.0
            do k = 1, N
                CM(i, j) = CM(i, j) + AM(i, k) * BM(k, j) / check
            end do
         end do
      end do
#else
! Loop 50
      do i = 1, N
         do  j = 1, N
            CM(i, j) = 0.0
            do k = 1, N
               if (i < j) then
                  CM(i, j) = CM(i, j) - AM(i, k) * BM(k, j) / check
               else
                  CM(i, j) = CM(i, j) + AM(i, k) * BM(k, j) / check
               endif
            end do
         end do
      end do
#endif

! Loop 60
#ifdef STRIDEACCESS
      do j = 1, N
         do i = 1, N
            sum = 0.0
            do k = 1, N
               sum = sum + CM(i, k) * AM (j, k)
            end do
            DM(i, j) = sum			
         end do
      end do

      do j = 1, N
        do i = 1, N
           CM(i, j) = DM(i, j)
        end do
      end do
#else
      do i = 1, N
         do j = 1, N
            sum = 0.0
            do k = 1, N
               sum = sum + CM(i, k) * AM (j, k)
            end do
            DM(i, j) = sum			
         end do
      end do

      do i = 1, N
        do j = 1, N
           CM(i, j) = DM(i, j)
        end do
      end do
#endif


! Loop 70
#ifdef STRIDEACCESS
	  do j = 1, N
		   do i = 1, N
			sum = 0.0
			do k = 1, N
			  sum = sum - CM(i, k) * BM (j, k)
			end do
			DM(i,j) = sum
		   end do
	  end do
#else
	  do i = 1, N
		   do j = 1, N
			sum = 0.0
			do k = 1, N
			  sum = sum - CM(i, k) * BM (j, k)
			end do
			DM(i,j) = sum
		   end do
	  end do
#endif


      HOLDA = abs(AM(1, 1))
      HOLDB = abs(BM(1, 1))
      do i = 1, N
        do j = 1, N
          HOLDA = max(HOLDA, abs(AM(i, j)))
          HOLDB = max(HOLDB, abs(BM(i, j)))
        end do
      end do

      TRACE3 = 0.0

! Loop 80
#ifdef TWO
        sub = HOLDA * HOLDB
#else
#endif
      do i = 1, N
#ifdef TWO
        TRACE3 = TRACE3 + (AM(IA(i),IA(i)) + BM(IA(i),IA(i)) - DM(IA(i),IA(i))) / sub
#else
        TRACE3 = TRACE3 + (AM(IA(i),IA(i)) + BM(IA(i),IA(i)) - DM(IA(i),IA(i))) / (HOLDA * HOLDB)
#endif
      end do

      call cpu_time(finish)
    print *, 'Final trace = ', trace3, ' and IDCHECK ', check
    print *, '-- RUNTIME -> ', finish-start, ' seconds'
    end


      double precision function trig (i,j)
      double precision x, y, z
#ifdef NOPRECISIONLOSS
      double precision pi 
#else
      real pi
#endif
      pi = acos(-1.0)
      x = dble(i) - dble(j)
      y = dble(i) + dble(j)
#ifdef ONE
      z = exp(sin(sqrt(x*x + y*y) * pi))
#else
      z = exp(sin(sqrt(x ** 2 + y ** 2) * pi))
#endif
      trig = x + y + log10(abs(1 + z + (x * y * z))) / (abs(x) + abs(y))
      return
      end

      subroutine idcheck(N,check,AV,BV,ID)

      double precision AV(*), BV(*), ID(N,*)
      double precision l2
      double precision check, check2
      double precision a, b, c, d, pastA, pastB, pastC, sub
      integer aux
	  
      real twoPi
      twoPi = 2.0 * acos(-1.0)

      pastA = 0
      pastB = 0
      pastC = 0
      check2 = 0
	  
#ifdef FOUR
#ifdef STRIDEACCESS
      do j = 1, N
        do i = 1, N
            ID(i,j) =  cos(check + i * twoPi / N) + 2.0 * sin(check + j * twoPi / N)
        end do
      end do
#else
      do i = 1, N
        do j = 1, N
            ID(i,j) =  cos(check + i * twoPi / N) + 2.0 * sin(check + j * twoPi / N)
        end do
      end do
#endif

      j = 1
      do i = 1, N
        if (((AV(i) < 0) .and. (BV(j) > 0 )) .or. ((AV(i) > 0) .and. (BV(j) < 0 ))) then
            ID(i,j) = -1.0
        else
            ID(i,j) = 1.0
        endif
        j = j + 1;
      end do
#else
#ifdef STRIDEACCESS
      do j = 1, N
        do i = 1, N
          if ( i == j ) then
             if ((AV(i) < 0) .and. (BV(j) < 0)) then
               ID(i,j) = 1.0
             elseif ((AV(i) < 0) .and. (BV(j) > 0 )) then
               ID(i,j) = -1.0
             elseif ((AV(i) > 0) .and. (BV(j) < 0 )) then
               ID(i,j) = -1.0
             else
               ID(i,j) = 1.0
             endif
          elseif (i /= j) then
            ID(i,j) =  cos(check + i * twoPi / N) + 2.0 * sin(check + j * twoPi / N)
          endif
        end do
      end do
#else
      do i = 1, N
        do j = 1, N
          if ( i == j ) then
             if ((AV(i) < 0) .and. (BV(j) < 0)) then
               ID(i,j) = 1.0
             elseif ((AV(i) < 0) .and. (BV(j) > 0 )) then
               ID(i,j) = -1.0
             elseif ((AV(i) > 0) .and. (BV(j) < 0 )) then
               ID(i,j) = -1.0
             else
               ID(i,j) = 1.0
             endif
          elseif (i /= j) then
            ID(i,j) =  cos(check + i * twoPi / N) + 2.0 * sin(check + j * twoPi / N)
          endif
        end do
      end do
#endif
#endif

      l2 = 0.0
#ifdef ONE
      do i = 1, N
        l2 = l2 + AV(i)*AV(i)
		!print *, "--- AV(i) ---", AV(i), "--- l2 ---", l2
      end do
#else
      do i = 1, N
        l2 = l2 + AV(i)**2
      end do
#endif

	  !print *, "--- l2 ---", l2, "--- sqrt(l2) ---", sqrt(l2)
      l2 = sqrt(l2)
      do i = 1, N
		!print *, "-->AV(i)", AV(i), "-->l2", l2 
        AV(i) = AV(i) / l2
      end do

      l2 = 0.0
      do i = 1, N
#ifdef ONE
      l2 = l2 + BV(i)*BV(i)
#else
      l2 = l2 + BV(i)**2
#endif
      end do

      l2 = sqrt(l2)
      do i = 1, N
        BV(i) = BV(i) / l2
      end do
	  
#ifdef FOUR
      a = 0.0
      b = 0.0
      c = 0.0
      d = 0.0  
! 0
          do i = 1, N
            do j = 1, N
#ifdef TWO
				sub = AV(i) * BV(j)
#else 
#endif 
              do k = 4 - int(mod(j + i, 4)), N, 4 
                pastA = a  
#ifdef TWO
                a = a + sub * ID(j, k)
#else
                a = a + AV(i) * BV(j) * ID(j, k)
#endif 
			end do
        end do
      end do
! 1
          do i = 1, N
            do j = 1, N
#ifdef TWO
            sub = AV(j) * BV(i)
#else
#endif 
          do k = 4 - int(mod(j + i - 1, 4)), N, 4
            pastB = b
#ifdef TWO
            b = b + sub * ID(k, j)
#else
            b = b + AV(j) * BV(i) * ID(k, j)
#endif
          end do
        end do
      end do

! 2
          do i = 1, N
            do j = 1, N
#ifdef TWO
            sub = AV(i) * BV(j)
#else
#endif
          do k = 4 - int(mod(j + i - 2, 4)), N, 4
            pastC = c
#ifdef TWO
            c = c - sub * ID(k, j)
#else
            c = c - AV(i) * BV(j) * ID(k, j)
#endif
          end do
        end do
      end do

! 3
          do i = 1, N
            do j = 1, N
#ifdef TWO
            sub = AV(j) * BV(i)
#else
#endif
          do k = 4 - int(mod(j + i + 1, 4)), N, 4
#ifdef TWO
            d = d - sub * ID(j, k)
#else
            d = d - AV(j) * BV(i) * ID(j, k)
#endif
          end do
        end do
      end do

      aux = 3 * N
      aux = int(mod(aux, 4))
      if (aux == 0) then
        check = sqrt(b*b + c*c) + a
        check2 = pastA + b + c + d
      else if (aux == 1) then
        check = sqrt(pastB*pastB + c*c) + a - b
        check2 = pastA + pastB + c + d
      else if (aux == 2) then
        check = sqrt(b*b + c*c)
        check2 = pastA + pastB + pastC + d
      else if (aux == 3) then
        check = sqrt(b*b + c*c)
        check2 = a + b + c + d
      end if
#else
      a = 0.0
      b = 0.0
      c = 0.0
      d = 0.0

      do i = 1, N
        do j = 1, N
          do k = 1, N
            aux = int(mod(i + j + k, 4))
            if (aux == 0) then
                a  = a +  AV(i) * BV(j) * ID(j, k)
                check = check + a
            else if (aux == 1) then
                b  = b +  AV(j) * BV(i) * ID(k, j)
                check = check - b
            else if (aux == 2) then
                c  = c -  AV(i) * BV(j) * ID(k, j)
#ifdef ONE
                check = sqrt(b*b + c*c)
#else
                check = sqrt(b**2 + c**2)
#endif
            else if (aux == 3) then
                d  = d -  AV(j) * BV(i) * ID(j, k)
                check2 = a + b + c + d
            else
            end if
          end do
        end do
      end do
#endif

		check = min(abs(check2), abs(check)) / max(abs(check2), abs(check))

      return
      end subroutine

    real (kind=8) function conrand(seed)
    !
    ! Function to generate a sequence of random numbers.
    ! Adapted from the  "Minimal Standard Method, real version 1 in Pascal"
    ! Park, S, Miller, K. "Random Number Generators: Good Ones are
    ! Hard to Find".  Communications of the ACM. vol 31, number 10,
    ! October 1988. pp. 1192-1201.
    !
    ! Fortran 2003 Version tested on 64 Bit Linux, gfortran compiler
    ! Andrew J. Pounds, Ph.D.
    ! Departments of Chemistry and Computer Science
    ! Mercer University
    ! Fall 2011
    !
        real (kind=8) :: seed
        real (kind=8) :: a, m
        real (kind=8) :: temp
        a = 16807.0D0
        m = 2147483647.0D0
        temp = a*seed
        seed = temp - m * int(temp/m)
        conrand = seed / m
    return
    end function conrand