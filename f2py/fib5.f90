      ! SUBROUTINE FIB5(T)

      ! INTEGER, parameter :: N = 100000000
      ! REAL*8 A(N), T
      ! integer(8) :: i, nsizemax, nsize, nloop, foo, cr, cm, tic, toc

      ! CALL SYSTEM_CLOCK(count_rate=cr, count_max=cm)
      ! CALL SYSTEM_CLOCK(tic)
      ! ! N = 100000000
      ! DO I=1,N
      !    IF (I.EQ.1) THEN
      !       A(I) = 0.0D0
      !    ELSEIF (I.EQ.2) THEN
      !       A(I) = 1.0D0
      !    ELSE 
      !       A(I) = A(I-1) + A(I-2)
      !    ENDIF
      ! ENDDO
      ! CALL SYSTEM_CLOCK(toc)

      ! T = (toc - tic) / dble(cr)

      ! ! print *, 'test'
      ! END


      SUBROUTINE FIB5(A)

      
      REAL(8), dimension(:), intent(inout) :: A
      INTEGER N
      integer(8) :: i, nsizemax, nsize, nloop, foo, cr, cm, tic, toc
      real(8) :: start, end

      ! CALL SYSTEM_CLOCK(count_rate=cr, count_max=cm)
      ! CALL SYSTEM_CLOCK(tic)
      ! CALL cpu_time(tic)

      N = size(A)
      DO I=1,N
         IF (I.EQ.1) THEN
            A(I) = 0.0D0
         ELSEIF (I.EQ.2) THEN
            A(I) = 1.0D0
         ELSEIF (I > 4) THEN
            A(I) = A(I-1) + A(I-2) + A(I-3) + A(I-4)
         ELSE 
            A(I) = A(I-1) + A(I-2)
         ENDIF
      ENDDO

      ! CALL SYSTEM_CLOCK(toc)

      ! T = (toc - tic) / dble(cr)

      ! write(*, A)
      END