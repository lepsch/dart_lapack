      SUBROUTINE DORHR_COL01( M, N, MB1, NB1, NB2, RESULT )
      IMPLICIT NONE

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int               M, N, MB1, NB1, NB2;
      // .. Return values ..
      double            RESULT(6);

*  =====================================================================

      // ..
      // .. Local allocatable arrays
      double          , ALLOCATABLE ::  A(:,:), AF(:,:), Q(:,:), R(:,:), RWORK(:), WORK( : ), T1(:,:), T2(:,:), DIAG(:), C(:,:), CF(:,:), D(:,:), DF(:,:);

      // .. Parameters ..
      double             ONE, ZERO;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               TESTZEROS;
      int                INFO, I, J, K, L, LWORK, NB1_UB, NB2_UB, NRB;
      double             ANORM, EPS, RESID, CNORM, DNORM;
      // ..
      // .. Local Arrays ..
      int                ISEED( 4 );
      double             WORKQUERY( 1 );
      // ..
      // .. External Functions ..
      double             DLAMCH, DLANGE, DLANSY;
      // EXTERNAL DLAMCH, DLANGE, DLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLACPY, DLARNV, DLASET, DLATSQR, DORHR_COL, DORGTSQR, DSCAL, DGEMM, DGEMQRT, DSYRK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CEILING, DBLE, MAX, MIN
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

      EPS = DLAMCH( 'Epsilon' )
      K = MIN( M, N )
      L = MAX( M, N, 1)

      // Dynamically allocate local arrays

      ALLOCATE ( A(M,N), AF(M,N), Q(L,L), R(M,L), RWORK(L), C(M,N), CF(M,N), D(N,M), DF(N,M) )

      // Put random numbers into A and copy to AF

      DO J = 1, N
         CALL DLARNV( 2, ISEED, M, A( 1, J ) )
      END DO
      if ( TESTZEROS ) {
         if ( M.GE.4 ) {
            DO J = 1, N
               CALL DLARNV( 2, ISEED, M/2, A( M/4, J ) )
            END DO
         }
      }
      CALL DLACPY( 'Full', M, N, A, M, AF, M )

      // Number of row blocks in DLATSQR

      NRB = MAX( 1, CEILING( DBLE( M - N ) / DBLE( MB1 - N ) ) )

      ALLOCATE ( T1( NB1, N * NRB ) )
      ALLOCATE ( T2( NB2, N ) )
      ALLOCATE ( DIAG( N ) )

      // Begin determine LWORK for the array WORK and allocate memory.

      // DLATSQR requires NB1 to be bounded by N.

      NB1_UB = MIN( NB1, N)

      // DGEMQRT requires NB2 to be bounded by N.

      NB2_UB = MIN( NB2, N)

      CALL DLATSQR( M, N, MB1, NB1_UB, AF, M, T1, NB1, WORKQUERY, -1, INFO )
      LWORK = INT( WORKQUERY( 1 ) )
      CALL DORGTSQR( M, N, MB1, NB1, AF, M, T1, NB1, WORKQUERY, -1, INFO )

      LWORK = MAX( LWORK, INT( WORKQUERY( 1 ) ) )

      // In DGEMQRT, WORK is N*NB2_UB if SIDE = 'L',
                 // or  M*NB2_UB if SIDE = 'R'.

      LWORK = MAX( LWORK, NB2_UB * N, NB2_UB * M )

      ALLOCATE ( WORK( LWORK ) )

      // End allocate memory for WORK.


      // Begin Householder reconstruction routines

      // Factor the matrix A in the array AF.

      SRNAMT = 'DLATSQR'
      CALL DLATSQR( M, N, MB1, NB1_UB, AF, M, T1, NB1, WORK, LWORK, INFO )

      // Copy the factor R into the array R.

      SRNAMT = 'DLACPY'
      CALL DLACPY( 'U', N, N, AF, M, R, M )

      // Reconstruct the orthogonal matrix Q.

      SRNAMT = 'DORGTSQR'
      CALL DORGTSQR( M, N, MB1, NB1, AF, M, T1, NB1, WORK, LWORK, INFO )

      // Perform the Householder reconstruction, the result is stored
     t // he arrays AF and T2.

      SRNAMT = 'DORHR_COL'
      CALL DORHR_COL( M, N, NB2, AF, M, T2, NB2, DIAG, INFO )

      // Compute the factor R_hr corresponding to the Householder
      // reconstructed Q_hr and place it in the upper triangle of AF to
      // match the Q storage format in DGEQRT. R_hr = R_tsqr * S,
     t // his means changing the sign of I-th row of the matrix R_tsqr
      // according to sign of of I-th diagonal element DIAG(I) of the
      // matrix S.

      SRNAMT = 'DLACPY'
      CALL DLACPY( 'U', N, N, R, M, AF, M )

      DO I = 1, N
         if ( DIAG( I ).EQ.-ONE ) {
            CALL DSCAL( N+1-I, -ONE, AF( I, I ), M )
         }
      END DO

      // End Householder reconstruction routines.


      // Generate the m-by-m matrix Q

      CALL DLASET( 'Full', M, M, ZERO, ONE, Q, M )

      SRNAMT = 'DGEMQRT'
      CALL DGEMQRT( 'L', 'N', M, M, K, NB2_UB, AF, M, T2, NB2, Q, M, WORK, INFO )

      // Copy R

      CALL DLASET( 'Full', M, N, ZERO, ZERO, R, M )

      CALL DLACPY( 'Upper', M, N, AF, M, R, M )

      // TEST 1
      // Compute |R - (Q**T)*A| / ( eps * m * |A| ) and store in RESULT(1)

      CALL DGEMM( 'T', 'N', M, N, M, -ONE, Q, M, A, M, ONE, R, M )

      ANORM = DLANGE( '1', M, N, A, M, RWORK )
      RESID = DLANGE( '1', M, N, R, M, RWORK )
      if ( ANORM.GT.ZERO ) {
         RESULT( 1 ) = RESID / ( EPS * MAX( 1, M ) * ANORM )
      } else {
         RESULT( 1 ) = ZERO
      }

      // TEST 2
      // Compute |I - (Q**T)*Q| / ( eps * m ) and store in RESULT(2)

      CALL DLASET( 'Full', M, M, ZERO, ONE, R, M )
      CALL DSYRK( 'U', 'T', M, M, -ONE, Q, M, ONE, R, M )
      RESID = DLANSY( '1', 'Upper', M, R, M, RWORK )
      RESULT( 2 ) = RESID / ( EPS * MAX( 1, M ) )

      // Generate random m-by-n matrix C

      DO J = 1, N
         CALL DLARNV( 2, ISEED, M, C( 1, J ) )
      END DO
      CNORM = DLANGE( '1', M, N, C, M, RWORK )
      CALL DLACPY( 'Full', M, N, C, M, CF, M )

      // Apply Q to C as Q*C = CF

      SRNAMT = 'DGEMQRT'
      CALL DGEMQRT( 'L', 'N', M, N, K, NB2_UB, AF, M, T2, NB2, CF, M, WORK, INFO )

      // TEST 3
      // Compute |CF - Q*C| / ( eps *  m * |C| )

      CALL DGEMM( 'N', 'N', M, N, M, -ONE, Q, M, C, M, ONE, CF, M )
      RESID = DLANGE( '1', M, N, CF, M, RWORK )
      if ( CNORM.GT.ZERO ) {
         RESULT( 3 ) = RESID / ( EPS * MAX( 1, M ) * CNORM )
      } else {
         RESULT( 3 ) = ZERO
      }

      // Copy C into CF again

      CALL DLACPY( 'Full', M, N, C, M, CF, M )

      // Apply Q to C as (Q**T)*C = CF

      SRNAMT = 'DGEMQRT'
      CALL DGEMQRT( 'L', 'T', M, N, K, NB2_UB, AF, M, T2, NB2, CF, M, WORK, INFO )

      // TEST 4
      // Compute |CF - (Q**T)*C| / ( eps * m * |C|)

      CALL DGEMM( 'T', 'N', M, N, M, -ONE, Q, M, C, M, ONE, CF, M )
      RESID = DLANGE( '1', M, N, CF, M, RWORK )
      if ( CNORM.GT.ZERO ) {
         RESULT( 4 ) = RESID / ( EPS * MAX( 1, M ) * CNORM )
      } else {
         RESULT( 4 ) = ZERO
      }

      // Generate random n-by-m matrix D and a copy DF

      DO J = 1, M
         CALL DLARNV( 2, ISEED, N, D( 1, J ) )
      END DO
      DNORM = DLANGE( '1', N, M, D, N, RWORK )
      CALL DLACPY( 'Full', N, M, D, N, DF, N )

      // Apply Q to D as D*Q = DF

      SRNAMT = 'DGEMQRT'
      CALL DGEMQRT( 'R', 'N', N, M, K, NB2_UB, AF, M, T2, NB2, DF, N, WORK, INFO )

      // TEST 5
      // Compute |DF - D*Q| / ( eps * m * |D| )

      CALL DGEMM( 'N', 'N', N, M, M, -ONE, D, N, Q, M, ONE, DF, N )
      RESID = DLANGE( '1', N, M, DF, N, RWORK )
      if ( DNORM.GT.ZERO ) {
         RESULT( 5 ) = RESID / ( EPS * MAX( 1, M ) * DNORM )
      } else {
         RESULT( 5 ) = ZERO
      }

      // Copy D into DF again

      CALL DLACPY( 'Full', N, M, D, N, DF, N )

      // Apply Q to D as D*QT = DF

      SRNAMT = 'DGEMQRT'
      CALL DGEMQRT( 'R', 'T', N, M, K, NB2_UB, AF, M, T2, NB2, DF, N, WORK, INFO )

      // TEST 6
      // Compute |DF - D*(Q**T)| / ( eps * m * |D| )

      CALL DGEMM( 'N', 'T', N, M, M, -ONE, D, N, Q, M, ONE, DF, N )
      RESID = DLANGE( '1', N, M, DF, N, RWORK )
      if ( DNORM.GT.ZERO ) {
         RESULT( 6 ) = RESID / ( EPS * MAX( 1, M ) * DNORM )
      } else {
         RESULT( 6 ) = ZERO
      }

      // Deallocate all arrays

      DEALLOCATE ( A, AF, Q, R, RWORK, WORK, T1, T2, DIAG, C, D, CF, DF )

      RETURN

      // End of DORHR_COL01

      }
