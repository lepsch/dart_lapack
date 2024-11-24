// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void slasd3(final int NL, final int NR, final int SQRE, final int K, final int D, final Matrix<double> Q_, final int LDQ, final int DSIGMA, final Matrix<double> U_, final int LDU, final Matrix<double> U2_, final int LDU2, final Matrix<double> VT_, final int LDVT, final Matrix<double> VT2_, final int LDVT2, final int IDXC, final int CTOT, final int Z, final Box<int> INFO,) {
  final Q = Q_.dim();
  final U = U_.dim();
  final U2 = U2_.dim();
  final VT = VT_.dim();
  final VT2 = VT2_.dim();

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, K, LDQ, LDU, LDU2, LDVT, LDVT2, NL, NR, SQRE;
      int                CTOT( * ), IDXC( * );
      double               D( * ), DSIGMA( * ), Q( LDQ, * ), U( LDU, * ), U2( LDU2, * ), VT( LDVT, * ), VT2( LDVT2, * ), Z( * );
      // ..

      double               ONE, ZERO, NEGONE;
      const              ONE = 1.0, ZERO = 0.0, NEGONE = -1.0 ;
      int                CTEMP, I, J, JC, KTEMP, M, N, NLP1, NLP2, NRP1;
      double               RHO, TEMP;
      // ..
      // .. External Functions ..
      //- REAL               SNRM2;
      // EXTERNAL SNRM2
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SGEMM, SLACPY, SLASCL, SLASD4, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, SIGN, SQRT

      // Test the input parameters.

      INFO = 0;

      if ( NL < 1 ) {
         INFO = -1;
      } else if ( NR < 1 ) {
         INFO = -2;
      } else if ( ( SQRE != 1 ) && ( SQRE != 0 ) ) {
         INFO = -3;
      }

      N = NL + NR + 1;
      M = N + SQRE;
      NLP1 = NL + 1;
      NLP2 = NL + 2;

      if ( ( K < 1 ) || ( K > N ) ) {
         INFO = -4;
      } else if ( LDQ < K ) {
         INFO = -7;
      } else if ( LDU < N ) {
         INFO = -10;
      } else if ( LDU2 < N ) {
         INFO = -12;
      } else if ( LDVT < M ) {
         INFO = -14;
      } else if ( LDVT2 < M ) {
         INFO = -16;
      }
      if ( INFO != 0 ) {
         xerbla('SLASD3', -INFO );
         return;
      }

      // Quick return if possible

      if ( K == 1 ) {
         D[1] = ( Z( 1 ) ).abs();
         scopy(M, VT2( 1, 1 ), LDVT2, VT( 1, 1 ), LDVT );
         if ( Z( 1 ) > ZERO ) {
            scopy(N, U2( 1, 1 ), 1, U( 1, 1 ), 1 );
         } else {
            for (I = 1; I <= N; I++) { // 10
               U[I][1] = -U2( I, 1 );
            } // 10
         }
         return;
      }

      // Keep a copy of Z.

      scopy(K, Z, 1, Q, 1 );

      // Normalize Z.

      RHO = SNRM2( K, Z, 1 );
      slascl('G', 0, 0, RHO, ONE, K, 1, Z, K, INFO );
      RHO = RHO*RHO;

      // Find the new singular values.

      for (J = 1; J <= K; J++) { // 30
         slasd4(K, J, DSIGMA, Z, U( 1, J ), RHO, D( J ), VT( 1, J ), INFO );

         // If the zero finder fails, report the convergence failure.

         if ( INFO != 0 ) {
            return;
         }
      } // 30

      // Compute updated Z.

      for (I = 1; I <= K; I++) { // 60
         Z[I] = U( I, K )*VT( I, K );
         for (J = 1; J <= I - 1; J++) { // 40
            Z[I] = Z( I )*( U( I, J )*VT( I, J ) / ( DSIGMA( I )-DSIGMA( J ) ) / ( DSIGMA( I )+DSIGMA( J ) ) );
         } // 40
         for (J = I; J <= K - 1; J++) { // 50
            Z[I] = Z( I )*( U( I, J )*VT( I, J ) / ( DSIGMA( I )-DSIGMA( J+1 ) ) / ( DSIGMA( I )+DSIGMA( J+1 ) ) );
         } // 50
         Z[I] = sign( sqrt( ( Z( I ) ).abs() ), Q( I, 1 ) );
      } // 60

      // Compute left singular vectors of the modified diagonal matrix,
      // and store related information for the right singular vectors.

      for (I = 1; I <= K; I++) { // 90
         VT[1][I] = Z( 1 ) / U( 1, I ) / VT( 1, I );
         U[1][I] = NEGONE;
         for (J = 2; J <= K; J++) { // 70
            VT[J][I] = Z( J ) / U( J, I ) / VT( J, I );
            U[J][I] = DSIGMA( J )*VT( J, I );
         } // 70
         TEMP = SNRM2( K, U( 1, I ), 1 );
         Q[1][I] = U( 1, I ) / TEMP;
         for (J = 2; J <= K; J++) { // 80
            JC = IDXC( J );
            Q[J][I] = U( JC, I ) / TEMP;
         } // 80
      } // 90

      // Update the left singular vector matrix.

      if ( K == 2 ) {
         sgemm('N', 'N', N, K, K, ONE, U2, LDU2, Q, LDQ, ZERO, U, LDU );
         GO TO 100;
      }
      if ( CTOT( 1 ) > 0 ) {
         sgemm('N', 'N', NL, K, CTOT( 1 ), ONE, U2( 1, 2 ), LDU2, Q( 2, 1 ), LDQ, ZERO, U( 1, 1 ), LDU );
         if ( CTOT( 3 ) > 0 ) {
            KTEMP = 2 + CTOT( 1 ) + CTOT( 2 );
            sgemm('N', 'N', NL, K, CTOT( 3 ), ONE, U2( 1, KTEMP ), LDU2, Q( KTEMP, 1 ), LDQ, ONE, U( 1, 1 ), LDU );
         }
      } else if ( CTOT( 3 ) > 0 ) {
         KTEMP = 2 + CTOT( 1 ) + CTOT( 2 );
         sgemm('N', 'N', NL, K, CTOT( 3 ), ONE, U2( 1, KTEMP ), LDU2, Q( KTEMP, 1 ), LDQ, ZERO, U( 1, 1 ), LDU );
      } else {
         slacpy('F', NL, K, U2, LDU2, U, LDU );
      }
      scopy(K, Q( 1, 1 ), LDQ, U( NLP1, 1 ), LDU );
      KTEMP = 2 + CTOT( 1 );
      CTEMP = CTOT( 2 ) + CTOT( 3 );
      sgemm('N', 'N', NR, K, CTEMP, ONE, U2( NLP2, KTEMP ), LDU2, Q( KTEMP, 1 ), LDQ, ZERO, U( NLP2, 1 ), LDU );

      // Generate the right singular vectors.

      } // 100
      for (I = 1; I <= K; I++) { // 120
         TEMP = SNRM2( K, VT( 1, I ), 1 );
         Q[I][1] = VT( 1, I ) / TEMP;
         for (J = 2; J <= K; J++) { // 110
            JC = IDXC( J );
            Q[I][J] = VT( JC, I ) / TEMP;
         } // 110
      } // 120

      // Update the right singular vector matrix.

      if ( K == 2 ) {
         sgemm('N', 'N', K, M, K, ONE, Q, LDQ, VT2, LDVT2, ZERO, VT, LDVT );
         return;
      }
      KTEMP = 1 + CTOT( 1 );
      sgemm('N', 'N', K, NLP1, KTEMP, ONE, Q( 1, 1 ), LDQ, VT2( 1, 1 ), LDVT2, ZERO, VT( 1, 1 ), LDVT );
      KTEMP = 2 + CTOT( 1 ) + CTOT( 2 );
      if (KTEMP <= LDVT2) sgemm( 'N', 'N', K, NLP1, CTOT( 3 ), ONE, Q( 1, KTEMP ), LDQ, VT2( KTEMP, 1 ), LDVT2, ONE, VT( 1, 1 ), LDVT );

      KTEMP = CTOT( 1 ) + 1;
      NRP1 = NR + SQRE;
      if ( KTEMP > 1 ) {
         for (I = 1; I <= K; I++) { // 130
            Q[I][KTEMP] = Q( I, 1 );
         } // 130
         for (I = NLP2; I <= M; I++) { // 140
            VT2[KTEMP][I] = VT2( 1, I );
         } // 140
      }
      CTEMP = 1 + CTOT( 2 ) + CTOT( 3 );
      sgemm('N', 'N', K, NRP1, CTEMP, ONE, Q( 1, KTEMP ), LDQ, VT2( KTEMP, NLP2 ), LDVT2, ZERO, VT( 1, NLP2 ), LDVT );

      }
