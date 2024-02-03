      SUBROUTINE ZUNHR_COL01( M, N, MB1, NB1, NB2, RESULT )
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
      COMPLEX*16      , ALLOCATABLE ::  A(:,:), AF(:,:), Q(:,:), R(:,:), WORK( : ), T1(:,:), T2(:,:), DIAG(:), C(:,:), CF(:,:), D(:,:), DF(:,:)
      double          , ALLOCATABLE :: RWORK(:);

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0D+0 ;
      COMPLEX*16         CONE, CZERO
      const              CONE = ( 1.0D+0, 0.0D+0 ), CZERO = ( 0.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               TESTZEROS;
      int                INFO, I, J, K, L, LWORK, NB1_UB, NB2_UB, NRB;
      double             ANORM, EPS, RESID, CNORM, DNORM;
      // ..
      // .. Local Arrays ..
      int                ISEED( 4 );
      COMPLEX*16         WORKQUERY( 1 )
      // ..
      // .. External Functions ..
      double             DLAMCH, ZLANGE, ZLANSY;
      // EXTERNAL DLAMCH, ZLANGE, ZLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLACPY, ZLARNV, ZLASET, ZLATSQR, ZUNHR_COL, ZUNGTSQR, ZSCAL, ZGEMM, ZGEMQRT, ZHERK
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
         CALL ZLARNV( 2, ISEED, M, A( 1, J ) )
      END DO
      if ( TESTZEROS ) {
         if ( M.GE.4 ) {
            DO J = 1, N
               CALL ZLARNV( 2, ISEED, M/2, A( M/4, J ) )
            END DO
         }
      }
      CALL ZLACPY( 'Full', M, N, A, M, AF, M )

      // Number of row blocks in ZLATSQR

      NRB = MAX( 1, CEILING( DBLE( M - N ) / DBLE( MB1 - N ) ) )

      ALLOCATE ( T1( NB1, N * NRB ) )
      ALLOCATE ( T2( NB2, N ) )
      ALLOCATE ( DIAG( N ) )

      // Begin determine LWORK for the array WORK and allocate memory.

      // ZLATSQR requires NB1 to be bounded by N.

      NB1_UB = MIN( NB1, N)

      // ZGEMQRT requires NB2 to be bounded by N.

      NB2_UB = MIN( NB2, N)

      CALL ZLATSQR( M, N, MB1, NB1_UB, AF, M, T1, NB1, WORKQUERY, -1, INFO )
      LWORK = INT( WORKQUERY( 1 ) )
      CALL ZUNGTSQR( M, N, MB1, NB1, AF, M, T1, NB1, WORKQUERY, -1, INFO )

      LWORK = MAX( LWORK, INT( WORKQUERY( 1 ) ) )

      // In ZGEMQRT, WORK is N*NB2_UB if SIDE = 'L',
                 // or  M*NB2_UB if SIDE = 'R'.

      LWORK = MAX( LWORK, NB2_UB * N, NB2_UB * M )

      ALLOCATE ( WORK( LWORK ) )

      // End allocate memory for WORK.


      // Begin Householder reconstruction routines

      // Factor the matrix A in the array AF.

      SRNAMT = 'ZLATSQR'
      CALL ZLATSQR( M, N, MB1, NB1_UB, AF, M, T1, NB1, WORK, LWORK, INFO )

      // Copy the factor R into the array R.

      SRNAMT = 'ZLACPY'
      CALL ZLACPY( 'U', N, N, AF, M, R, M )

      // Reconstruct the orthogonal matrix Q.

      SRNAMT = 'ZUNGTSQR'
      CALL ZUNGTSQR( M, N, MB1, NB1, AF, M, T1, NB1, WORK, LWORK, INFO )

      // Perform the Householder reconstruction, the result is stored
     t // he arrays AF and T2.

      SRNAMT = 'ZUNHR_COL'
      CALL ZUNHR_COL( M, N, NB2, AF, M, T2, NB2, DIAG, INFO )

      // Compute the factor R_hr corresponding to the Householder
      // reconstructed Q_hr and place it in the upper triangle of AF to
      // match the Q storage format in ZGEQRT. R_hr = R_tsqr * S,
     t // his means changing the sign of I-th row of the matrix R_tsqr
      // according to sign of of I-th diagonal element DIAG(I) of the
      // matrix S.

      SRNAMT = 'ZLACPY'
      CALL ZLACPY( 'U', N, N, R, M, AF, M )

      DO I = 1, N
         if ( DIAG( I ).EQ.-CONE ) {
            CALL ZSCAL( N+1-I, -CONE, AF( I, I ), M )
         }
      END DO

      // End Householder reconstruction routines.


      // Generate the m-by-m matrix Q

      CALL ZLASET( 'Full', M, M, CZERO, CONE, Q, M )

      SRNAMT = 'ZGEMQRT'
      CALL ZGEMQRT( 'L', 'N', M, M, K, NB2_UB, AF, M, T2, NB2, Q, M, WORK, INFO )

      // Copy R

      CALL ZLASET( 'Full', M, N, CZERO, CZERO, R, M )

      CALL ZLACPY( 'Upper', M, N, AF, M, R, M )

      // TEST 1
      // Compute |R - (Q**H)*A| / ( eps * m * |A| ) and store in RESULT(1)

      CALL ZGEMM( 'C', 'N', M, N, M, -CONE, Q, M, A, M, CONE, R, M )

      ANORM = ZLANGE( '1', M, N, A, M, RWORK )
      RESID = ZLANGE( '1', M, N, R, M, RWORK )
      if ( ANORM.GT.ZERO ) {
         RESULT( 1 ) = RESID / ( EPS * MAX( 1, M ) * ANORM )
      } else {
         RESULT( 1 ) = ZERO
      }

      // TEST 2
      // Compute |I - (Q**H)*Q| / ( eps * m ) and store in RESULT(2)

      CALL ZLASET( 'Full', M, M, CZERO, CONE, R, M )
      CALL ZHERK( 'U', 'C', M, M, REAL(-CONE), Q, M, REAL(CONE), R, M )
      RESID = ZLANSY( '1', 'Upper', M, R, M, RWORK )
      RESULT( 2 ) = RESID / ( EPS * MAX( 1, M ) )

      // Generate random m-by-n matrix C

      DO J = 1, N
         CALL ZLARNV( 2, ISEED, M, C( 1, J ) )
      END DO
      CNORM = ZLANGE( '1', M, N, C, M, RWORK )
      CALL ZLACPY( 'Full', M, N, C, M, CF, M )

      // Apply Q to C as Q*C = CF

      SRNAMT = 'ZGEMQRT'
      CALL ZGEMQRT( 'L', 'N', M, N, K, NB2_UB, AF, M, T2, NB2, CF, M, WORK, INFO )

      // TEST 3
      // Compute |CF - Q*C| / ( eps *  m * |C| )

      CALL ZGEMM( 'N', 'N', M, N, M, -CONE, Q, M, C, M, CONE, CF, M )
      RESID = ZLANGE( '1', M, N, CF, M, RWORK )
      if ( CNORM.GT.ZERO ) {
         RESULT( 3 ) = RESID / ( EPS * MAX( 1, M ) * CNORM )
      } else {
         RESULT( 3 ) = ZERO
      }

      // Copy C into CF again

      CALL ZLACPY( 'Full', M, N, C, M, CF, M )

      // Apply Q to C as (Q**H)*C = CF

      SRNAMT = 'ZGEMQRT'
      CALL ZGEMQRT( 'L', 'C', M, N, K, NB2_UB, AF, M, T2, NB2, CF, M, WORK, INFO )

      // TEST 4
      // Compute |CF - (Q**H)*C| / ( eps * m * |C|)

      CALL ZGEMM( 'C', 'N', M, N, M, -CONE, Q, M, C, M, CONE, CF, M )
      RESID = ZLANGE( '1', M, N, CF, M, RWORK )
      if ( CNORM.GT.ZERO ) {
         RESULT( 4 ) = RESID / ( EPS * MAX( 1, M ) * CNORM )
      } else {
         RESULT( 4 ) = ZERO
      }

      // Generate random n-by-m matrix D and a copy DF

      DO J = 1, M
         CALL ZLARNV( 2, ISEED, N, D( 1, J ) )
      END DO
      DNORM = ZLANGE( '1', N, M, D, N, RWORK )
      CALL ZLACPY( 'Full', N, M, D, N, DF, N )

      // Apply Q to D as D*Q = DF

      SRNAMT = 'ZGEMQRT'
      CALL ZGEMQRT( 'R', 'N', N, M, K, NB2_UB, AF, M, T2, NB2, DF, N, WORK, INFO )

      // TEST 5
      // Compute |DF - D*Q| / ( eps * m * |D| )

      CALL ZGEMM( 'N', 'N', N, M, M, -CONE, D, N, Q, M, CONE, DF, N )
      RESID = ZLANGE( '1', N, M, DF, N, RWORK )
      if ( DNORM.GT.ZERO ) {
         RESULT( 5 ) = RESID / ( EPS * MAX( 1, M ) * DNORM )
      } else {
         RESULT( 5 ) = ZERO
      }

      // Copy D into DF again

      CALL ZLACPY( 'Full', N, M, D, N, DF, N )

      // Apply Q to D as D*QT = DF

      SRNAMT = 'ZGEMQRT'
      CALL ZGEMQRT( 'R', 'C', N, M, K, NB2_UB, AF, M, T2, NB2, DF, N, WORK, INFO )

      // TEST 6
      // Compute |DF - D*(Q**H)| / ( eps * m * |D| )

      CALL ZGEMM( 'N', 'C', N, M, M, -CONE, D, N, Q, M, CONE, DF, N )
      RESID = ZLANGE( '1', N, M, DF, N, RWORK )
      if ( DNORM.GT.ZERO ) {
         RESULT( 6 ) = RESID / ( EPS * MAX( 1, M ) * DNORM )
      } else {
         RESULT( 6 ) = ZERO
      }

      // Deallocate all arrays

      DEALLOCATE ( A, AF, Q, R, RWORK, WORK, T1, T2, DIAG, C, D, CF, DF )

      RETURN

      // End of ZUNHR_COL01

      }
