      SUBROUTINE SORHR_COL01( M, N, MB1, NB1, NB2, RESULT )
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
      int                INFO, I, J, K, L, LWORK, NB1_UB, NB2_UB, NRB;
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
      // EXTERNAL SLACPY, SLARNV, SLASET, SLATSQR, SORHR_COL, SORGTSQR, SSCAL, SGEMM, SGEMQRT, SSYRK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CEILING, REAL, MAX, MIN
      // ..
      // .. Scalars in Common ..
      String   (LEN=32)  SRNAMT;
      // ..
      // .. Common blocks ..
      // COMMON / SRMNAMC / SRNAMT
      // ..
      // .. Data statements ..
      DATA ISEED / 1988, 1989, 1990, 1991 /

      // TEST MATRICES WITH HALF OF MATRIX BEING ZEROS

      TESTZEROS = false;

      EPS = SLAMCH( 'Epsilon' )
      K = MIN( M, N )
      L = MAX( M, N, 1)

      // Dynamically allocate local arrays

      ALLOCATE ( A(M,N), AF(M,N), Q(L,L), R(M,L), RWORK(L), C(M,N), CF(M,N), D(N,M), DF(N,M) )

      // Put random numbers into A and copy to AF

      for (J = 1; J <= N; J++) {
         slarnv(2, ISEED, M, A( 1, J ) );
      }
      if ( TESTZEROS ) {
         if ( M >= 4 ) {
            for (J = 1; J <= N; J++) {
               slarnv(2, ISEED, M/2, A( M/4, J ) );
            }
         }
      }
      slacpy('Full', M, N, A, M, AF, M );

      // Number of row blocks in SLATSQR

      NRB = MAX( 1, CEILING( REAL( M - N ) / REAL( MB1 - N ) ) )

      ALLOCATE ( T1( NB1, N * NRB ) )
      ALLOCATE ( T2( NB2, N ) )
      ALLOCATE ( DIAG( N ) )

      // Begin determine LWORK for the array WORK and allocate memory.

      // SLATSQR requires NB1 to be bounded by N.

      NB1_UB = MIN( NB1, N)

      // SGEMQRT requires NB2 to be bounded by N.

      NB2_UB = MIN( NB2, N)

      slatsqr(M, N, MB1, NB1_UB, AF, M, T1, NB1, WORKQUERY, -1, INFO );
      LWORK = INT( WORKQUERY( 1 ) )
      sorgtsqr(M, N, MB1, NB1, AF, M, T1, NB1, WORKQUERY, -1, INFO );

      LWORK = MAX( LWORK, INT( WORKQUERY( 1 ) ) )

      // In SGEMQRT, WORK is N*NB2_UB if SIDE = 'L',
                 // or  M*NB2_UB if SIDE = 'R'.

      LWORK = MAX( LWORK, NB2_UB * N, NB2_UB * M )

      ALLOCATE ( WORK( LWORK ) )

      // End allocate memory for WORK.


      // Begin Householder reconstruction routines

      // Factor the matrix A in the array AF.

      SRNAMT = 'SLATSQR'
      slatsqr(M, N, MB1, NB1_UB, AF, M, T1, NB1, WORK, LWORK, INFO );

      // Copy the factor R into the array R.

      SRNAMT = 'SLACPY'
      slacpy('U', N, N, AF, M, R, M );

      // Reconstruct the orthogonal matrix Q.

      SRNAMT = 'SORGTSQR'
      sorgtsqr(M, N, MB1, NB1, AF, M, T1, NB1, WORK, LWORK, INFO );

      // Perform the Householder reconstruction, the result is stored
      // the arrays AF and T2.

      SRNAMT = 'SORHR_COL'
      sorhr_col(M, N, NB2, AF, M, T2, NB2, DIAG, INFO );

      // Compute the factor R_hr corresponding to the Householder
      // reconstructed Q_hr and place it in the upper triangle of AF to
      // match the Q storage format in SGEQRT. R_hr = R_tsqr * S,
      // this means changing the sign of I-th row of the matrix R_tsqr
      // according to sign of of I-th diagonal element DIAG(I) of the
      // matrix S.

      SRNAMT = 'SLACPY'
      slacpy('U', N, N, R, M, AF, M );

      for (I = 1; I <= N; I++) {
         if ( DIAG( I ) == -ONE ) {
            sscal(N+1-I, -ONE, AF( I, I ), M );
         }
      }

      // End Householder reconstruction routines.


      // Generate the m-by-m matrix Q

      slaset('Full', M, M, ZERO, ONE, Q, M );

      SRNAMT = 'SGEMQRT'
      sgemqrt('L', 'N', M, M, K, NB2_UB, AF, M, T2, NB2, Q, M, WORK, INFO );

      // Copy R

      slaset('Full', M, N, ZERO, ZERO, R, M );

      slacpy('Upper', M, N, AF, M, R, M );

      // TEST 1
      // Compute |R - (Q**T)*A| / ( eps * m * |A| ) and store in RESULT(1)

      sgemm('T', 'N', M, N, M, -ONE, Q, M, A, M, ONE, R, M );

      ANORM = SLANGE( '1', M, N, A, M, RWORK )
      RESID = SLANGE( '1', M, N, R, M, RWORK )
      if ( ANORM > ZERO ) {
         RESULT( 1 ) = RESID / ( EPS * MAX( 1, M ) * ANORM )
      } else {
         RESULT( 1 ) = ZERO
      }

      // TEST 2
      // Compute |I - (Q**T)*Q| / ( eps * m ) and store in RESULT(2)

      slaset('Full', M, M, ZERO, ONE, R, M );
      ssyrk('U', 'T', M, M, -ONE, Q, M, ONE, R, M );
      RESID = SLANSY( '1', 'Upper', M, R, M, RWORK )
      RESULT( 2 ) = RESID / ( EPS * MAX( 1, M ) )

      // Generate random m-by-n matrix C

      for (J = 1; J <= N; J++) {
         slarnv(2, ISEED, M, C( 1, J ) );
      }
      CNORM = SLANGE( '1', M, N, C, M, RWORK )
      slacpy('Full', M, N, C, M, CF, M );

      // Apply Q to C as Q*C = CF

      SRNAMT = 'SGEMQRT'
      sgemqrt('L', 'N', M, N, K, NB2_UB, AF, M, T2, NB2, CF, M, WORK, INFO );

      // TEST 3
      // Compute |CF - Q*C| / ( eps *  m * |C| )

      sgemm('N', 'N', M, N, M, -ONE, Q, M, C, M, ONE, CF, M );
      RESID = SLANGE( '1', M, N, CF, M, RWORK )
      if ( CNORM > ZERO ) {
         RESULT( 3 ) = RESID / ( EPS * MAX( 1, M ) * CNORM )
      } else {
         RESULT( 3 ) = ZERO
      }

      // Copy C into CF again

      slacpy('Full', M, N, C, M, CF, M );

      // Apply Q to C as (Q**T)*C = CF

      SRNAMT = 'SGEMQRT'
      sgemqrt('L', 'T', M, N, K, NB2_UB, AF, M, T2, NB2, CF, M, WORK, INFO );

      // TEST 4
      // Compute |CF - (Q**T)*C| / ( eps * m * |C|)

      sgemm('T', 'N', M, N, M, -ONE, Q, M, C, M, ONE, CF, M );
      RESID = SLANGE( '1', M, N, CF, M, RWORK )
      if ( CNORM > ZERO ) {
         RESULT( 4 ) = RESID / ( EPS * MAX( 1, M ) * CNORM )
      } else {
         RESULT( 4 ) = ZERO
      }

      // Generate random n-by-m matrix D and a copy DF

      for (J = 1; J <= M; J++) {
         slarnv(2, ISEED, N, D( 1, J ) );
      }
      DNORM = SLANGE( '1', N, M, D, N, RWORK )
      slacpy('Full', N, M, D, N, DF, N );

      // Apply Q to D as D*Q = DF

      SRNAMT = 'SGEMQRT'
      sgemqrt('R', 'N', N, M, K, NB2_UB, AF, M, T2, NB2, DF, N, WORK, INFO );

      // TEST 5
      // Compute |DF - D*Q| / ( eps * m * |D| )

      sgemm('N', 'N', N, M, M, -ONE, D, N, Q, M, ONE, DF, N );
      RESID = SLANGE( '1', N, M, DF, N, RWORK )
      if ( DNORM > ZERO ) {
         RESULT( 5 ) = RESID / ( EPS * MAX( 1, M ) * DNORM )
      } else {
         RESULT( 5 ) = ZERO
      }

      // Copy D into DF again

      slacpy('Full', N, M, D, N, DF, N );

      // Apply Q to D as D*QT = DF

      SRNAMT = 'SGEMQRT'
      sgemqrt('R', 'T', N, M, K, NB2_UB, AF, M, T2, NB2, DF, N, WORK, INFO );

      // TEST 6
      // Compute |DF - D*(Q**T)| / ( eps * m * |D| )

      sgemm('N', 'T', N, M, M, -ONE, D, N, Q, M, ONE, DF, N );
      RESID = SLANGE( '1', N, M, DF, N, RWORK )
      if ( DNORM > ZERO ) {
         RESULT( 6 ) = RESID / ( EPS * MAX( 1, M ) * DNORM )
      } else {
         RESULT( 6 ) = ZERO
      }

      // Deallocate all arrays

      DEALLOCATE ( A, AF, Q, R, RWORK, WORK, T1, T2, DIAG, C, D, CF, DF )

      RETURN

      // End of SORHR_COL01

      }
