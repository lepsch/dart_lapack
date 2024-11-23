      SUBROUTINE SORHR_COL02( M, N, MB1, NB1, NB2, RESULT )
      IMPLICIT NONE

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int               M, N, MB1, NB1, NB2;
      // .. Return values ..
      REAL              RESULT(6)

*  =====================================================================

      // ..
      // .. Local allocatable arrays
      REAL            , ALLOCATABLE ::  A(:,:), AF(:,:), Q(:,:), R(:,:), RWORK(:), WORK( : ), T1(:,:), T2(:,:), DIAG(:), C(:,:), CF(:,:), D(:,:), DF(:,:)

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               TESTZEROS;
      int                INFO, J, K, L, LWORK, NB2_UB, NRB;
      REAL               ANORM, EPS, RESID, CNORM, DNORM
      // ..
      // .. Local Arrays ..
      int                ISEED( 4 );
      REAL               WORKQUERY( 1 )
      // ..
      // .. External Functions ..
      REAL               SLAMCH, SLANGE, SLANSY
      // EXTERNAL SLAMCH, SLANGE, SLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLACPY, SLARNV, SLASET, SGETSQRHRT, SSCAL, SGEMM, SGEMQRT, SSYRK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CEILING, REAL, MAX, MIN
      // ..
      // .. Scalars in Common ..
      String   (LEN=32)  SRNAMT;
      // ..
      // .. Common blocks ..
      COMMON             / SRMNAMC / SRNAMT
      // ..
      // .. Data statements ..
      DATA ISEED / 1988, 1989, 1990, 1991 /

      // TEST MATRICES WITH HALF OF MATRIX BEING ZEROS

      TESTZEROS = .FALSE.

      EPS = SLAMCH( 'Epsilon' )
      K = MIN( M, N )
      L = MAX( M, N, 1)

      // Dynamically allocate local arrays

      ALLOCATE ( A(M,N), AF(M,N), Q(L,L), R(M,L), RWORK(L), C(M,N), CF(M,N), D(N,M), DF(N,M) )

      // Put random numbers into A and copy to AF

      DO J = 1, N
         CALL SLARNV( 2, ISEED, M, A( 1, J ) )
      END DO
      IF( TESTZEROS ) THEN
         IF( M.GE.4 ) THEN
            DO J = 1, N
               CALL SLARNV( 2, ISEED, M/2, A( M/4, J ) )
            END DO
         END IF
      END IF
      CALL SLACPY( 'Full', M, N, A, M, AF, M )

      // Number of row blocks in SLATSQR

      NRB = MAX( 1, CEILING( REAL( M - N ) / REAL( MB1 - N ) ) )

      ALLOCATE ( T1( NB1, N * NRB ) )
      ALLOCATE ( T2( NB2, N ) )
      ALLOCATE ( DIAG( N ) )

      // Begin determine LWORK for the array WORK and allocate memory.

      // SGEMQRT requires NB2 to be bounded by N.

      NB2_UB = MIN( NB2, N)

      CALL SGETSQRHRT( M, N, MB1, NB1, NB2, AF, M, T2, NB2, WORKQUERY, -1, INFO )

      LWORK = INT( WORKQUERY( 1 ) )

      // In SGEMQRT, WORK is N*NB2_UB if SIDE = 'L',
                 // or  M*NB2_UB if SIDE = 'R'.

      LWORK = MAX( LWORK, NB2_UB * N, NB2_UB * M )

      ALLOCATE ( WORK( LWORK ) )

      // End allocate memory for WORK.


      // Begin Householder reconstruction routines

      // Factor the matrix A in the array AF.

      SRNAMT = 'SGETSQRHRT'
      CALL SGETSQRHRT( M, N, MB1, NB1, NB2, AF, M, T2, NB2, WORK, LWORK, INFO )

      // End Householder reconstruction routines.


      // Generate the m-by-m matrix Q

      CALL SLASET( 'Full', M, M, ZERO, ONE, Q, M )

      SRNAMT = 'SGEMQRT'
      CALL SGEMQRT( 'L', 'N', M, M, K, NB2_UB, AF, M, T2, NB2, Q, M, WORK, INFO )

      // Copy R

      CALL SLASET( 'Full', M, N, ZERO, ZERO, R, M )

      CALL SLACPY( 'Upper', M, N, AF, M, R, M )

      // TEST 1
      // Compute |R - (Q**T)*A| / ( eps * m * |A| ) and store in RESULT(1)

      CALL SGEMM( 'T', 'N', M, N, M, -ONE, Q, M, A, M, ONE, R, M )

      ANORM = SLANGE( '1', M, N, A, M, RWORK )
      RESID = SLANGE( '1', M, N, R, M, RWORK )
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = RESID / ( EPS * MAX( 1, M ) * ANORM )
      } else {
         RESULT( 1 ) = ZERO
      END IF

      // TEST 2
      // Compute |I - (Q**T)*Q| / ( eps * m ) and store in RESULT(2)

      CALL SLASET( 'Full', M, M, ZERO, ONE, R, M )
      CALL SSYRK( 'U', 'T', M, M, -ONE, Q, M, ONE, R, M )
      RESID = SLANSY( '1', 'Upper', M, R, M, RWORK )
      RESULT( 2 ) = RESID / ( EPS * MAX( 1, M ) )

      // Generate random m-by-n matrix C

      DO J = 1, N
         CALL SLARNV( 2, ISEED, M, C( 1, J ) )
      END DO
      CNORM = SLANGE( '1', M, N, C, M, RWORK )
      CALL SLACPY( 'Full', M, N, C, M, CF, M )

      // Apply Q to C as Q*C = CF

      SRNAMT = 'SGEMQRT'
      CALL SGEMQRT( 'L', 'N', M, N, K, NB2_UB, AF, M, T2, NB2, CF, M, WORK, INFO )

      // TEST 3
      // Compute |CF - Q*C| / ( eps *  m * |C| )

      CALL SGEMM( 'N', 'N', M, N, M, -ONE, Q, M, C, M, ONE, CF, M )
      RESID = SLANGE( '1', M, N, CF, M, RWORK )
      IF( CNORM.GT.ZERO ) THEN
         RESULT( 3 ) = RESID / ( EPS * MAX( 1, M ) * CNORM )
      } else {
         RESULT( 3 ) = ZERO
      END IF

      // Copy C into CF again

      CALL SLACPY( 'Full', M, N, C, M, CF, M )

      // Apply Q to C as (Q**T)*C = CF

      SRNAMT = 'SGEMQRT'
      CALL SGEMQRT( 'L', 'T', M, N, K, NB2_UB, AF, M, T2, NB2, CF, M, WORK, INFO )

      // TEST 4
      // Compute |CF - (Q**T)*C| / ( eps * m * |C|)

      CALL SGEMM( 'T', 'N', M, N, M, -ONE, Q, M, C, M, ONE, CF, M )
      RESID = SLANGE( '1', M, N, CF, M, RWORK )
      IF( CNORM.GT.ZERO ) THEN
         RESULT( 4 ) = RESID / ( EPS * MAX( 1, M ) * CNORM )
      } else {
         RESULT( 4 ) = ZERO
      END IF

      // Generate random n-by-m matrix D and a copy DF

      DO J = 1, M
         CALL SLARNV( 2, ISEED, N, D( 1, J ) )
      END DO
      DNORM = SLANGE( '1', N, M, D, N, RWORK )
      CALL SLACPY( 'Full', N, M, D, N, DF, N )

      // Apply Q to D as D*Q = DF

      SRNAMT = 'SGEMQRT'
      CALL SGEMQRT( 'R', 'N', N, M, K, NB2_UB, AF, M, T2, NB2, DF, N, WORK, INFO )

      // TEST 5
      // Compute |DF - D*Q| / ( eps * m * |D| )

      CALL SGEMM( 'N', 'N', N, M, M, -ONE, D, N, Q, M, ONE, DF, N )
      RESID = SLANGE( '1', N, M, DF, N, RWORK )
      IF( DNORM.GT.ZERO ) THEN
         RESULT( 5 ) = RESID / ( EPS * MAX( 1, M ) * DNORM )
      } else {
         RESULT( 5 ) = ZERO
      END IF

      // Copy D into DF again

      CALL SLACPY( 'Full', N, M, D, N, DF, N )

      // Apply Q to D as D*QT = DF

      SRNAMT = 'SGEMQRT'
      CALL SGEMQRT( 'R', 'T', N, M, K, NB2_UB, AF, M, T2, NB2, DF, N, WORK, INFO )

      // TEST 6
      // Compute |DF - D*(Q**T)| / ( eps * m * |D| )

      CALL SGEMM( 'N', 'T', N, M, M, -ONE, D, N, Q, M, ONE, DF, N )
      RESID = SLANGE( '1', N, M, DF, N, RWORK )
      IF( DNORM.GT.ZERO ) THEN
         RESULT( 6 ) = RESID / ( EPS * MAX( 1, M ) * DNORM )
      } else {
         RESULT( 6 ) = ZERO
      END IF

      // Deallocate all arrays

      DEALLOCATE ( A, AF, Q, R, RWORK, WORK, T1, T2, DIAG, C, D, CF, DF )

      RETURN

      // End of SORHR_COL02

      }
