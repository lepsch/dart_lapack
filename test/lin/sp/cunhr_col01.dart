// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void cunhr_col01(final int M, final int N, final int MB1, final int NB1, final int NB2, final int RESULT,) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int               M, N, MB1, NB1, NB2;
      // .. Return values ..
      double              RESULT(6);

// =====================================================================

      // ..
      // .. Local allocatable arrays
      Complex         , ALLOCATABLE ::  A(:,:), AF(:,:), Q(:,:), R(:,:), WORK( : ), T1(:,:), T2(:,:), DIAG(:), C(:,:), CF(:,:), D(:,:), DF(:,:);
      double            , ALLOCATABLE :: RWORK(:);

      // .. Parameters ..
      double               ZERO;
      const              ZERO = 0.0 ;
      Complex            CONE, CZERO;
      const              CONE = ( 1.0, 0.0 ), CZERO = ( 0.0, 0.0 ) ;
      bool               TESTZEROS;
      int                INFO, I, J, K, L, LWORK, NB1_UB, NB2_UB, NRB;
      double               ANORM, EPS, RESID, CNORM, DNORM;
      int                ISEED( 4 );
      Complex            WORKQUERY( 1 );
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH, CLANGE, CLANSY;
      // EXTERNAL SLAMCH, CLANGE, CLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLACPY, CLARNV, CLASET, CLATSQR, CUNHR_COL, CUNGTSQR, CSCAL, CGEMM, CGEMQRT, CHERK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CEILING, REAL, MAX, MIN
      // ..
      // .. Scalars in Common ..
      String   (LEN=32) srnamc.SRNAMT;
      // ..
      // .. Common blocks ..
      // COMMON / SRMNAMC /srnamc.SRNAMT
      // ..
      // .. Data statements ..
      const ISEED = [ 1988, 1989, 1990, 1991 ];

      // TEST MATRICES WITH HALF OF MATRIX BEING ZEROS

      TESTZEROS = false;

      EPS = SLAMCH( 'Epsilon' );
      K = min( M, N );
      L = max( M, N, 1);

      // Dynamically allocate local arrays

      ALLOCATE ( A(M,N), AF(M,N), Q(L,L), R(M,L), RWORK(L), C(M,N), CF(M,N), D(N,M), DF(N,M) );

      // Put random numbers into A and copy to AF

      for (J = 1; J <= N; J++) {
         clarnv(2, ISEED, M, A( 1, J ) );
      }
      if ( TESTZEROS ) {
         if ( M >= 4 ) {
            for (J = 1; J <= N; J++) {
               clarnv(2, ISEED, M/2, A( M/4, J ) );
            }
         }
      }
      clacpy('Full', M, N, A, M, AF, M );

      // Number of row blocks in CLATSQR

      NRB = max( 1, CEILING( double( M - N ) / REAL( MB1 - N ) ) );

      ALLOCATE ( T1( NB1, N * NRB ) );
      ALLOCATE ( T2( NB2, N ) );
      ALLOCATE ( DIAG( N ) );

      // Begin determine LWORK for the array WORK and allocate memory.

      // CLATSQR requires NB1 to be bounded by N.

      NB1_UB = min( NB1, N);

      // CGEMQRT requires NB2 to be bounded by N.

      NB2_UB = min( NB2, N);

      clatsqr(M, N, MB1, NB1_UB, AF, M, T1, NB1, WORKQUERY, -1, INFO );
      LWORK = INT( WORKQUERY( 1 ) );
      cungtsqr(M, N, MB1, NB1, AF, M, T1, NB1, WORKQUERY, -1, INFO );

      LWORK = max( LWORK, INT( WORKQUERY( 1 ) ) );

      // In CGEMQRT, WORK is N*NB2_UB if SIDE = 'L',
      //            or  M*NB2_UB if SIDE = 'R'.

      LWORK = max( LWORK, NB2_UB * N, NB2_UB * M );

      ALLOCATE ( WORK( LWORK ) );

      // End allocate memory for WORK.


      // Begin Householder reconstruction routines

      // Factor the matrix A in the array AF.

     srnamc.SRNAMT = 'CLATSQR';
      clatsqr(M, N, MB1, NB1_UB, AF, M, T1, NB1, WORK, LWORK, INFO );

      // Copy the factor R into the array R.

     srnamc.SRNAMT = 'CLACPY';
      clacpy('U', N, N, AF, M, R, M );

      // Reconstruct the orthogonal matrix Q.

     srnamc.SRNAMT = 'CUNGTSQR';
      cungtsqr(M, N, MB1, NB1, AF, M, T1, NB1, WORK, LWORK, INFO );

      // Perform the Householder reconstruction, the result is stored
      // the arrays AF and T2.

     srnamc.SRNAMT = 'CUNHR_COL';
      cunhr_col(M, N, NB2, AF, M, T2, NB2, DIAG, INFO );

      // Compute the factor R_hr corresponding to the Householder
      // reconstructed Q_hr and place it in the upper triangle of AF to
      // match the Q storage format in CGEQRT. R_hr = R_tsqr * S,
      // this means changing the sign of I-th row of the matrix R_tsqr
      // according to sign of of I-th diagonal element DIAG(I) of the
      // matrix S.

     srnamc.SRNAMT = 'CLACPY';
      clacpy('U', N, N, R, M, AF, M );

      for (I = 1; I <= N; I++) {
         if ( DIAG( I ) == -CONE ) {
            cscal(N+1-I, -CONE, AF( I, I ), M );
         }
      }

      // End Householder reconstruction routines.


      // Generate the m-by-m matrix Q

      claset('Full', M, M, CZERO, CONE, Q, M );

     srnamc.SRNAMT = 'CGEMQRT';
      cgemqrt('L', 'N', M, M, K, NB2_UB, AF, M, T2, NB2, Q, M, WORK, INFO );

      // Copy R

      claset('Full', M, N, CZERO, CZERO, R, M );

      clacpy('Upper', M, N, AF, M, R, M );

      // TEST 1
      // Compute |R - (Q**H)*A| / ( eps * m * |A| ) and store in RESULT(1)

      cgemm('C', 'N', M, N, M, -CONE, Q, M, A, M, CONE, R, M );

      ANORM = CLANGE( '1', M, N, A, M, RWORK );
      RESID = CLANGE( '1', M, N, R, M, RWORK );
      if ( ANORM > ZERO ) {
         RESULT[1] = RESID / ( EPS * max( 1, M ) * ANORM );
      } else {
         RESULT[1] = ZERO;
      }

      // TEST 2
      // Compute |I - (Q**H)*Q| / ( eps * m ) and store in RESULT(2)

      claset('Full', M, M, CZERO, CONE, R, M );
      cherk('U', 'C', M, M, REAL(-CONE), Q, M, double(CONE), R, M );
      RESID = CLANSY( '1', 'Upper', M, R, M, RWORK );
      RESULT[2] = RESID / ( EPS * max( 1, M ) );

      // Generate random m-by-n matrix C

      for (J = 1; J <= N; J++) {
         clarnv(2, ISEED, M, C( 1, J ) );
      }
      CNORM = CLANGE( '1', M, N, C, M, RWORK );
      clacpy('Full', M, N, C, M, CF, M );

      // Apply Q to C as Q*C = CF

     srnamc.SRNAMT = 'CGEMQRT';
      cgemqrt('L', 'N', M, N, K, NB2_UB, AF, M, T2, NB2, CF, M, WORK, INFO );

      // TEST 3
      // Compute |CF - Q*C| / ( eps *  m * |C| )

      cgemm('N', 'N', M, N, M, -CONE, Q, M, C, M, CONE, CF, M );
      RESID = CLANGE( '1', M, N, CF, M, RWORK );
      if ( CNORM > ZERO ) {
         RESULT[3] = RESID / ( EPS * max( 1, M ) * CNORM );
      } else {
         RESULT[3] = ZERO;
      }

      // Copy C into CF again

      clacpy('Full', M, N, C, M, CF, M );

      // Apply Q to C as (Q**H)*C = CF

     srnamc.SRNAMT = 'CGEMQRT';
      cgemqrt('L', 'C', M, N, K, NB2_UB, AF, M, T2, NB2, CF, M, WORK, INFO );

      // TEST 4
      // Compute |CF - (Q**H)*C| / ( eps * m * |C|)

      cgemm('C', 'N', M, N, M, -CONE, Q, M, C, M, CONE, CF, M );
      RESID = CLANGE( '1', M, N, CF, M, RWORK );
      if ( CNORM > ZERO ) {
         RESULT[4] = RESID / ( EPS * max( 1, M ) * CNORM );
      } else {
         RESULT[4] = ZERO;
      }

      // Generate random n-by-m matrix D and a copy DF

      for (J = 1; J <= M; J++) {
         clarnv(2, ISEED, N, D( 1, J ) );
      }
      DNORM = CLANGE( '1', N, M, D, N, RWORK );
      clacpy('Full', N, M, D, N, DF, N );

      // Apply Q to D as D*Q = DF

     srnamc.SRNAMT = 'CGEMQRT';
      cgemqrt('R', 'N', N, M, K, NB2_UB, AF, M, T2, NB2, DF, N, WORK, INFO );

      // TEST 5
      // Compute |DF - D*Q| / ( eps * m * |D| )

      cgemm('N', 'N', N, M, M, -CONE, D, N, Q, M, CONE, DF, N );
      RESID = CLANGE( '1', N, M, DF, N, RWORK );
      if ( DNORM > ZERO ) {
         RESULT[5] = RESID / ( EPS * max( 1, M ) * DNORM );
      } else {
         RESULT[5] = ZERO;
      }

      // Copy D into DF again

      clacpy('Full', N, M, D, N, DF, N );

      // Apply Q to D as D*QT = DF

     srnamc.SRNAMT = 'CGEMQRT';
      cgemqrt('R', 'C', N, M, K, NB2_UB, AF, M, T2, NB2, DF, N, WORK, INFO );

      // TEST 6
      // Compute |DF - D*(Q**H)| / ( eps * m * |D| )

      cgemm('N', 'C', N, M, M, -CONE, D, N, Q, M, CONE, DF, N );
      RESID = CLANGE( '1', N, M, DF, N, RWORK );
      if ( DNORM > ZERO ) {
         RESULT[6] = RESID / ( EPS * max( 1, M ) * DNORM );
      } else {
         RESULT[6] = ZERO;
      }

      // Deallocate all arrays

      DEALLOCATE ( A, AF, Q, R, RWORK, WORK, T1, T2, DIAG, C, D, CF, DF );

      }
