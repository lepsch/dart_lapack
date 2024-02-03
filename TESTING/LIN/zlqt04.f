      SUBROUTINE ZLQT04(M,N,NB,RESULT)
      IMPLICIT NONE

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     M, N, NB;
      // .. Return values ..
      double           RESULT(6);

*  =====================================================================

      // ..
      // .. Local allocatable arrays
      COMPLEX*16, ALLOCATABLE :: AF(:,:), Q(:,:), L(:,:), WORK( : ), T(:,:), CF(:,:), DF(:,:), A(:,:), C(:,:), D(:,:)
      double          , ALLOCATABLE :: RWORK(:);

      // .. Parameters ..
      double           ZERO;
      COMPLEX*16 ONE, CZERO
      const    ZERO = 0.0;
      const    ONE = (1.0,0.0), CZERO=(0.0,0.0) ;
      // ..
      // .. Local Scalars ..
      int     INFO, J, K, LL, LWORK, LDT;
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
      LL = MAX(M,N)
      LWORK = MAX(2,LL)*MAX(2,LL)*NB

      // Dynamically allocate local arrays

      ALLOCATE ( A(M,N), AF(M,N), Q(N,N), L(LL,N), RWORK(LL), WORK(LWORK), T(NB,N), C(M,N), CF(M,N), D(N,M), DF(N,M) )

      // Put random numbers into A and copy to AF

      LDT=NB
      DO J=1,N
         CALL ZLARNV( 2, ISEED, M, A( 1, J ) )
      END DO
      CALL ZLACPY( 'Full', M, N, A, M, AF, M )

      // Factor the matrix A in the array AF.

      CALL ZGELQT( M, N, NB, AF, M, T, LDT, WORK, INFO )

      // Generate the n-by-n matrix Q

      CALL ZLASET( 'Full', N, N, CZERO, ONE, Q, N )
      CALL ZGEMLQT( 'R', 'N', N, N, K, NB, AF, M, T, LDT, Q, N, WORK, INFO )

      // Copy L

      CALL ZLASET( 'Full', LL, N, CZERO, CZERO, L, LL )
      CALL ZLACPY( 'Lower', M, N, AF, M, L, LL )

      // Compute |L - A*Q'| / |A| and store in RESULT(1)

      CALL ZGEMM( 'N', 'C', M, N, N, -ONE, A, M, Q, N, ONE, L, LL )
      ANORM = ZLANGE( '1', M, N, A, M, RWORK )
      RESID = ZLANGE( '1', M, N, L, LL, RWORK )
      if ( ANORM.GT.ZERO ) {
         RESULT( 1 ) = RESID / (EPS*MAX(1,M)*ANORM)
      } else {
         RESULT( 1 ) = ZERO
      }

      // Compute |I - Q'*Q| and store in RESULT(2)

      CALL ZLASET( 'Full', N, N, CZERO, ONE, L, LL )
      CALL ZHERK( 'U', 'C', N, N, DREAL(-ONE), Q, N, DREAL(ONE), L, LL)
      RESID = ZLANSY( '1', 'Upper', N, L, LL, RWORK )
      RESULT( 2 ) = RESID / (EPS*MAX(1,N))

      // Generate random m-by-n matrix C and a copy CF

      DO J=1,M
         CALL ZLARNV( 2, ISEED, N, D( 1, J ) )
      END DO
      DNORM = ZLANGE( '1', N, M, D, N, RWORK)
      CALL ZLACPY( 'Full', N, M, D, N, DF, N )

      // Apply Q to C as Q*C

      CALL ZGEMLQT( 'L', 'N', N, M, K, NB, AF, M, T, NB, DF, N, WORK, INFO)

      // Compute |Q*D - Q*D| / |D|

      CALL ZGEMM( 'N', 'N', N, M, N, -ONE, Q, N, D, N, ONE, DF, N )
      RESID = ZLANGE( '1', N, M, DF, N, RWORK )
      if ( DNORM.GT.ZERO ) {
         RESULT( 3 ) = RESID / (EPS*MAX(1,M)*DNORM)
      } else {
         RESULT( 3 ) = ZERO
      }

      // Copy D into DF again

      CALL ZLACPY( 'Full', N, M, D, N, DF, N )

      // Apply Q to D as QT*D

      CALL ZGEMLQT( 'L', 'C', N, M, K, NB, AF, M, T, NB, DF, N, WORK, INFO)

      // Compute |QT*D - QT*D| / |D|

      CALL ZGEMM( 'C', 'N', N, M, N, -ONE, Q, N, D, N, ONE, DF, N )
      RESID = ZLANGE( '1', N, M, DF, N, RWORK )
      if ( DNORM.GT.ZERO ) {
         RESULT( 4 ) = RESID / (EPS*MAX(1,M)*DNORM)
      } else {
         RESULT( 4 ) = ZERO
      }

      // Generate random n-by-m matrix D and a copy DF

      DO J=1,N
         CALL ZLARNV( 2, ISEED, M, C( 1, J ) )
      END DO
      CNORM = ZLANGE( '1', M, N, C, M, RWORK)
      CALL ZLACPY( 'Full', M, N, C, M, CF, M )

      // Apply Q to C as C*Q

      CALL ZGEMLQT( 'R', 'N', M, N, K, NB, AF, M, T, NB, CF, M, WORK, INFO)

      // Compute |C*Q - C*Q| / |C|

      CALL ZGEMM( 'N', 'N', M, N, N, -ONE, C, M, Q, N, ONE, CF, M )
      RESID = ZLANGE( '1', N, M, DF, N, RWORK )
      if ( CNORM.GT.ZERO ) {
         RESULT( 5 ) = RESID / (EPS*MAX(1,M)*DNORM)
      } else {
         RESULT( 5 ) = ZERO
      }

      // Copy C into CF again

      CALL ZLACPY( 'Full', M, N, C, M, CF, M )

      // Apply Q to D as D*QT

      CALL ZGEMLQT( 'R', 'C', M, N, K, NB, AF, M, T, NB, CF, M, WORK, INFO)

      // Compute |C*QT - C*QT| / |C|

      CALL ZGEMM( 'N', 'C', M, N, N, -ONE, C, M, Q, N, ONE, CF, M )
      RESID = ZLANGE( '1', M, N, CF, M, RWORK )
      if ( CNORM.GT.ZERO ) {
         RESULT( 6 ) = RESID / (EPS*MAX(1,M)*DNORM)
      } else {
         RESULT( 6 ) = ZERO
      }

      // Deallocate all arrays

      DEALLOCATE ( A, AF, Q, L, RWORK, WORK, T, C, D, CF, DF)

      RETURN
      }
