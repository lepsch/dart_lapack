      SUBROUTINE ZQRT04(M,N,NB,RESULT)
      IMPLICIT NONE

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     M, N, NB, LDT;
      // .. Return values ..
      double           RESULT(6);

*  =====================================================================

      // ..
      // .. Local allocatable arrays
      COMPLEX*16, ALLOCATABLE :: AF(:,:), Q(:,:), R(:,:), WORK( : ), T(:,:), CF(:,:), DF(:,:), A(:,:), C(:,:), D(:,:)
      double          , ALLOCATABLE :: RWORK(:);

      // .. Parameters ..
      double           ZERO;
      COMPLEX*16 ONE, CZERO
      const    ZERO = 0.0, ONE = (1.0,0.0), CZERO=(0.0,0.0) ;
      // ..
      // .. Local Scalars ..
      int     INFO, J, K, L, LWORK;
      double             ANORM, EPS, RESID, CNORM, DNORM;
      // ..
      // .. Local Arrays ..
      int                ISEED( 4 );
      // ..
      // .. External Functions ..
      double           DLAMCH;
      double           ZLANGE, ZLANSY;
      bool     LSAME;
      // EXTERNAL DLAMCH, ZLANGE, ZLANSY, LSAME
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Data statements ..
      DATA ISEED / 1988, 1989, 1990, 1991 /

      EPS = DLAMCH( 'Epsilon' )
      K = MIN(M,N)
      L = MAX(M,N)
      LWORK = MAX(2,L)*MAX(2,L)*NB

      // Dynamically allocate local arrays

      ALLOCATE ( A(M,N), AF(M,N), Q(M,M), R(M,L), RWORK(L), WORK(LWORK), T(NB,N), C(M,N), CF(M,N), D(N,M), DF(N,M) )

      // Put random numbers into A and copy to AF

      LDT=NB
      DO J=1,N
         CALL ZLARNV( 2, ISEED, M, A( 1, J ) )
      END DO
      CALL ZLACPY( 'Full', M, N, A, M, AF, M )

      // Factor the matrix A in the array AF.

      CALL ZGEQRT( M, N, NB, AF, M, T, LDT, WORK, INFO )

      // Generate the m-by-m matrix Q

      CALL ZLASET( 'Full', M, M, CZERO, ONE, Q, M )
      CALL ZGEMQRT( 'R', 'N', M, M, K, NB, AF, M, T, LDT, Q, M, WORK, INFO )

      // Copy R

      CALL ZLASET( 'Full', M, N, CZERO, CZERO, R, M )
      CALL ZLACPY( 'Upper', M, N, AF, M, R, M )

      // Compute |R - Q'*A| / |A| and store in RESULT(1)

      CALL ZGEMM( 'C', 'N', M, N, M, -ONE, Q, M, A, M, ONE, R, M )
      ANORM = ZLANGE( '1', M, N, A, M, RWORK )
      RESID = ZLANGE( '1', M, N, R, M, RWORK )
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = RESID / (EPS*MAX(1,M)*ANORM)
      ELSE
         RESULT( 1 ) = ZERO
      END IF

      // Compute |I - Q'*Q| and store in RESULT(2)

      CALL ZLASET( 'Full', M, M, CZERO, ONE, R, M )
      CALL ZHERK( 'U', 'C', M, M, DREAL(-ONE), Q, M, DREAL(ONE), R, M )
      RESID = ZLANSY( '1', 'Upper', M, R, M, RWORK )
      RESULT( 2 ) = RESID / (EPS*MAX(1,M))

      // Generate random m-by-n matrix C and a copy CF

      DO J=1,N
         CALL ZLARNV( 2, ISEED, M, C( 1, J ) )
      END DO
      CNORM = ZLANGE( '1', M, N, C, M, RWORK)
      CALL ZLACPY( 'Full', M, N, C, M, CF, M )

      // Apply Q to C as Q*C

      CALL ZGEMQRT( 'L', 'N', M, N, K, NB, AF, M, T, NB, CF, M, WORK, INFO)

      // Compute |Q*C - Q*C| / |C|

      CALL ZGEMM( 'N', 'N', M, N, M, -ONE, Q, M, C, M, ONE, CF, M )
      RESID = ZLANGE( '1', M, N, CF, M, RWORK )
      IF( CNORM.GT.ZERO ) THEN
         RESULT( 3 ) = RESID / (EPS*MAX(1,M)*CNORM)
      ELSE
         RESULT( 3 ) = ZERO
      END IF

      // Copy C into CF again

      CALL ZLACPY( 'Full', M, N, C, M, CF, M )

      // Apply Q to C as QT*C

      CALL ZGEMQRT( 'L', 'C', M, N, K, NB, AF, M, T, NB, CF, M, WORK, INFO)

      // Compute |QT*C - QT*C| / |C|

      CALL ZGEMM( 'C', 'N', M, N, M, -ONE, Q, M, C, M, ONE, CF, M )
      RESID = ZLANGE( '1', M, N, CF, M, RWORK )
      IF( CNORM.GT.ZERO ) THEN
         RESULT( 4 ) = RESID / (EPS*MAX(1,M)*CNORM)
      ELSE
         RESULT( 4 ) = ZERO
      END IF

      // Generate random n-by-m matrix D and a copy DF

      DO J=1,M
         CALL ZLARNV( 2, ISEED, N, D( 1, J ) )
      END DO
      DNORM = ZLANGE( '1', N, M, D, N, RWORK)
      CALL ZLACPY( 'Full', N, M, D, N, DF, N )

      // Apply Q to D as D*Q

      CALL ZGEMQRT( 'R', 'N', N, M, K, NB, AF, M, T, NB, DF, N, WORK, INFO)

      // Compute |D*Q - D*Q| / |D|

      CALL ZGEMM( 'N', 'N', N, M, M, -ONE, D, N, Q, M, ONE, DF, N )
      RESID = ZLANGE( '1', N, M, DF, N, RWORK )
      IF( CNORM.GT.ZERO ) THEN
         RESULT( 5 ) = RESID / (EPS*MAX(1,M)*DNORM)
      ELSE
         RESULT( 5 ) = ZERO
      END IF

      // Copy D into DF again

      CALL ZLACPY( 'Full', N, M, D, N, DF, N )

      // Apply Q to D as D*QT

      CALL ZGEMQRT( 'R', 'C', N, M, K, NB, AF, M, T, NB, DF, N, WORK, INFO)

      // Compute |D*QT - D*QT| / |D|

      CALL ZGEMM( 'N', 'C', N, M, M, -ONE, D, N, Q, M, ONE, DF, N )
      RESID = ZLANGE( '1', N, M, DF, N, RWORK )
      IF( CNORM.GT.ZERO ) THEN
         RESULT( 6 ) = RESID / (EPS*MAX(1,M)*DNORM)
      ELSE
         RESULT( 6 ) = ZERO
      END IF

      // Deallocate all arrays

      DEALLOCATE ( A, AF, Q, R, RWORK, WORK, T, C, D, CF, DF)

      RETURN
      }
