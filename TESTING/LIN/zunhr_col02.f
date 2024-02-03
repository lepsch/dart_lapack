      SUBROUTINE ZUNHR_COL02( M, N, MB1, NB1, NB2, RESULT )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER           M, N, MB1, NB1, NB2
*     .. Return values ..
      DOUBLE PRECISION  RESULT(6)
*
*  =====================================================================
*
*     ..
*     .. Local allocatable arrays
      COMPLEX*16      , ALLOCATABLE ::  A(:,:), AF(:,:), Q(:,:), R(:,:), WORK( : ), T1(:,:), T2(:,:), DIAG(:), C(:,:), CF(:,:), D(:,:), DF(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: RWORK(:)
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
      COMPLEX*16         CONE, CZERO
      PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ), CZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            TESTZEROS
      INTEGER            INFO, J, K, L, LWORK, NB2_UB, NRB
      DOUBLE PRECISION   ANORM, EPS, RESID, CNORM, DNORM
*     ..
*     .. Local Arrays ..
      INTEGER            ISEED( 4 )
      COMPLEX*16         WORKQUERY( 1 )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, ZLANGE, ZLANSY
      EXTERNAL           DLAMCH, ZLANGE, ZLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZLACPY, ZLARNV, ZLASET, ZGETSQRHRT, ZSCAL, ZGEMM, ZGEMQRT, ZHERK
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CEILING, DBLE, MAX, MIN
*     ..
*     .. Scalars in Common ..
      CHARACTER(LEN=32)  SRNAMT
*     ..
*     .. Common blocks ..
      COMMON             / SRMNAMC / SRNAMT
*     ..
*     .. Data statements ..
      DATA ISEED / 1988, 1989, 1990, 1991 /
*
*     TEST MATRICES WITH HALF OF MATRIX BEING ZEROS
*
      TESTZEROS = .FALSE.
*
      EPS = DLAMCH( 'Epsilon' )
      K = MIN( M, N )
      L = MAX( M, N, 1)
*
*     Dynamically allocate local arrays
*
      ALLOCATE ( A(M,N), AF(M,N), Q(L,L), R(M,L), RWORK(L), C(M,N), CF(M,N), D(N,M), DF(N,M) )
*
*     Put random numbers into A and copy to AF
*
      DO J = 1, N
         CALL ZLARNV( 2, ISEED, M, A( 1, J ) )
      END DO
      IF( TESTZEROS ) THEN
         IF( M.GE.4 ) THEN
            DO J = 1, N
               CALL ZLARNV( 2, ISEED, M/2, A( M/4, J ) )
            END DO
         END IF
      END IF
      CALL ZLACPY( 'Full', M, N, A, M, AF, M )
*
*     Number of row blocks in ZLATSQR
*
      NRB = MAX( 1, CEILING( DBLE( M - N ) / DBLE( MB1 - N ) ) )
*
      ALLOCATE ( T1( NB1, N * NRB ) )
      ALLOCATE ( T2( NB2, N ) )
      ALLOCATE ( DIAG( N ) )
*
*     Begin determine LWORK for the array WORK and allocate memory.
*
*     ZGEMQRT requires NB2 to be bounded by N.
*
      NB2_UB = MIN( NB2, N)
*
*
      CALL ZGETSQRHRT( M, N, MB1, NB1, NB2, AF, M, T2, NB2, WORKQUERY, -1, INFO )
*
      LWORK = INT( WORKQUERY( 1 ) )
*
*     In ZGEMQRT, WORK is N*NB2_UB if SIDE = 'L',
*                or  M*NB2_UB if SIDE = 'R'.
*
      LWORK = MAX( LWORK, NB2_UB * N, NB2_UB * M )
*
      ALLOCATE ( WORK( LWORK ) )
*
*     End allocate memory for WORK.
*
*
*     Begin Householder reconstruction routines
*
*     Factor the matrix A in the array AF.
*
      SRNAMT = 'ZGETSQRHRT'
      CALL ZGETSQRHRT( M, N, MB1, NB1, NB2, AF, M, T2, NB2, WORK, LWORK, INFO )
*
*     End Householder reconstruction routines.
*
*
*     Generate the m-by-m matrix Q
*
      CALL ZLASET( 'Full', M, M, CZERO, CONE, Q, M )
*
      SRNAMT = 'ZGEMQRT'
      CALL ZGEMQRT( 'L', 'N', M, M, K, NB2_UB, AF, M, T2, NB2, Q, M, WORK, INFO )
*
*     Copy R
*
      CALL ZLASET( 'Full', M, N, CZERO, CZERO, R, M )
*
      CALL ZLACPY( 'Upper', M, N, AF, M, R, M )
*
*     TEST 1
*     Compute |R - (Q**T)*A| / ( eps * m * |A| ) and store in RESULT(1)
*
      CALL ZGEMM( 'C', 'N', M, N, M, -CONE, Q, M, A, M, CONE, R, M )
*
      ANORM = ZLANGE( '1', M, N, A, M, RWORK )
      RESID = ZLANGE( '1', M, N, R, M, RWORK )
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = RESID / ( EPS * MAX( 1, M ) * ANORM )
      ELSE
         RESULT( 1 ) = ZERO
      END IF
*
*     TEST 2
*     Compute |I - (Q**T)*Q| / ( eps * m ) and store in RESULT(2)
*
      CALL ZLASET( 'Full', M, M, CZERO, CONE, R, M )
      CALL ZHERK( 'U', 'C', M, M, REAL(-CONE), Q, M, REAL(CONE), R, M )
      RESID = ZLANSY( '1', 'Upper', M, R, M, RWORK )
      RESULT( 2 ) = RESID / ( EPS * MAX( 1, M ) )
*
*     Generate random m-by-n matrix C
*
      DO J = 1, N
         CALL ZLARNV( 2, ISEED, M, C( 1, J ) )
      END DO
      CNORM = ZLANGE( '1', M, N, C, M, RWORK )
      CALL ZLACPY( 'Full', M, N, C, M, CF, M )
*
*     Apply Q to C as Q*C = CF
*
      SRNAMT = 'ZGEMQRT'
      CALL ZGEMQRT( 'L', 'N', M, N, K, NB2_UB, AF, M, T2, NB2, CF, M, WORK, INFO )
*
*     TEST 3
*     Compute |CF - Q*C| / ( eps *  m * |C| )
*
      CALL ZGEMM( 'N', 'N', M, N, M, -CONE, Q, M, C, M, CONE, CF, M )
      RESID = ZLANGE( '1', M, N, CF, M, RWORK )
      IF( CNORM.GT.ZERO ) THEN
         RESULT( 3 ) = RESID / ( EPS * MAX( 1, M ) * CNORM )
      ELSE
         RESULT( 3 ) = ZERO
      END IF
*
*     Copy C into CF again
*
      CALL ZLACPY( 'Full', M, N, C, M, CF, M )
*
*     Apply Q to C as (Q**T)*C = CF
*
      SRNAMT = 'ZGEMQRT'
      CALL ZGEMQRT( 'L', 'C', M, N, K, NB2_UB, AF, M, T2, NB2, CF, M, WORK, INFO )
*
*     TEST 4
*     Compute |CF - (Q**T)*C| / ( eps * m * |C|)
*
      CALL ZGEMM( 'C', 'N', M, N, M, -CONE, Q, M, C, M, CONE, CF, M )
      RESID = ZLANGE( '1', M, N, CF, M, RWORK )
      IF( CNORM.GT.ZERO ) THEN
         RESULT( 4 ) = RESID / ( EPS * MAX( 1, M ) * CNORM )
      ELSE
         RESULT( 4 ) = ZERO
      END IF
*
*     Generate random n-by-m matrix D and a copy DF
*
      DO J = 1, M
         CALL ZLARNV( 2, ISEED, N, D( 1, J ) )
      END DO
      DNORM = ZLANGE( '1', N, M, D, N, RWORK )
      CALL ZLACPY( 'Full', N, M, D, N, DF, N )
*
*     Apply Q to D as D*Q = DF
*
      SRNAMT = 'ZGEMQRT'
      CALL ZGEMQRT( 'R', 'N', N, M, K, NB2_UB, AF, M, T2, NB2, DF, N, WORK, INFO )
*
*     TEST 5
*     Compute |DF - D*Q| / ( eps * m * |D| )
*
      CALL ZGEMM( 'N', 'N', N, M, M, -CONE, D, N, Q, M, CONE, DF, N )
      RESID = ZLANGE( '1', N, M, DF, N, RWORK )
      IF( DNORM.GT.ZERO ) THEN
         RESULT( 5 ) = RESID / ( EPS * MAX( 1, M ) * DNORM )
      ELSE
         RESULT( 5 ) = ZERO
      END IF
*
*     Copy D into DF again
*
      CALL ZLACPY( 'Full', N, M, D, N, DF, N )
*
*     Apply Q to D as D*QT = DF
*
      SRNAMT = 'ZGEMQRT'
      CALL ZGEMQRT( 'R', 'C', N, M, K, NB2_UB, AF, M, T2, NB2, DF, N, WORK, INFO )
*
*     TEST 6
*     Compute |DF - D*(Q**T)| / ( eps * m * |D| )
*
      CALL ZGEMM( 'N', 'C', N, M, M, -CONE, D, N, Q, M, CONE, DF, N )
      RESID = ZLANGE( '1', N, M, DF, N, RWORK )
      IF( DNORM.GT.ZERO ) THEN
         RESULT( 6 ) = RESID / ( EPS * MAX( 1, M ) * DNORM )
      ELSE
         RESULT( 6 ) = ZERO
      END IF
*
*     Deallocate all arrays
*
      DEALLOCATE ( A, AF, Q, R, RWORK, WORK, T1, T2, DIAG, C, D, CF, DF )
*
      RETURN
*
*     End of ZUNHR_COL02
*
      END
