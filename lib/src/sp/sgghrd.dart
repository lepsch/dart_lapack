// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void sgghrd(final int COMPQ, final int COMPZ, final int N, final int ILO, final int IHI, final Matrix<double> A_, final int LDA, final Matrix<double> B_, final int LDB, final Matrix<double> Q_, final int LDQ, final Matrix<double> Z_, final int LDZ, final Box<int> INFO,) {
  final A = A_.dim();
  final B = B_.dim();
  final Q = Q_.dim();
  final Z = Z_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             COMPQ, COMPZ;
      int                IHI, ILO, INFO, LDA, LDB, LDQ, LDZ, N;
      double               A( LDA, * ), B( LDB, * ), Q( LDQ, * ), Z( LDZ, * );
      // ..

      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      bool               ILQ, ILZ;
      int                ICOMPQ, ICOMPZ, JCOL, JROW;
      double               C, S, TEMP;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLARTG, SLASET, SROT, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      // Decode COMPQ

      if ( lsame( COMPQ, 'N' ) ) {
         ILQ = false;
         ICOMPQ = 1;
      } else if ( lsame( COMPQ, 'V' ) ) {
         ILQ = true;
         ICOMPQ = 2;
      } else if ( lsame( COMPQ, 'I' ) ) {
         ILQ = true;
         ICOMPQ = 3;
      } else {
         ICOMPQ = 0;
      }

      // Decode COMPZ

      if ( lsame( COMPZ, 'N' ) ) {
         ILZ = false;
         ICOMPZ = 1;
      } else if ( lsame( COMPZ, 'V' ) ) {
         ILZ = true;
         ICOMPZ = 2;
      } else if ( lsame( COMPZ, 'I' ) ) {
         ILZ = true;
         ICOMPZ = 3;
      } else {
         ICOMPZ = 0;
      }

      // Test the input parameters.

      INFO = 0;
      if ( ICOMPQ <= 0 ) {
         INFO = -1;
      } else if ( ICOMPZ <= 0 ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( ILO < 1 ) {
         INFO = -4;
      } else if ( IHI > N || IHI < ILO-1 ) {
         INFO = -5;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -7;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -9;
      } else if ( ( ILQ && LDQ < N ) || LDQ < 1 ) {
         INFO = -11;
      } else if ( ( ILZ && LDZ < N ) || LDZ < 1 ) {
         INFO = -13;
      }
      if ( INFO != 0 ) {
         xerbla('SGGHRD', -INFO );
         return;
      }

      // Initialize Q and Z if desired.

      if (ICOMPQ == 3) slaset( 'Full', N, N, ZERO, ONE, Q, LDQ );
      IF( ICOMPZ == 3 ) slaset( 'Full', N, N, ZERO, ONE, Z, LDZ );

      // Quick return if possible

      if (N <= 1) return;

      // Zero out lower triangle of B

      for (JCOL = 1; JCOL <= N - 1; JCOL++) { // 20
         for (JROW = JCOL + 1; JROW <= N; JROW++) { // 10
            B[JROW][JCOL] = ZERO;
         } // 10
      } // 20

      // Reduce A and B

      for (JCOL = ILO; JCOL <= IHI - 2; JCOL++) { // 40

         for (JROW = IHI; JROW >= JCOL + 2; JROW--) { // 30

            // Step 1: rotate rows JROW-1, JROW to kill A(JROW,JCOL)

            TEMP = A( JROW-1, JCOL );
            slartg(TEMP, A( JROW, JCOL ), C, S, A( JROW-1, JCOL ) );
            A[JROW][JCOL] = ZERO;
            srot(N-JCOL, A( JROW-1, JCOL+1 ), LDA, A( JROW, JCOL+1 ), LDA, C, S );
            srot(N+2-JROW, B( JROW-1, JROW-1 ), LDB, B( JROW, JROW-1 ), LDB, C, S )             IF( ILQ ) CALL SROT( N, Q( 1, JROW-1 ), 1, Q( 1, JROW ), 1, C, S );

            // Step 2: rotate columns JROW, JROW-1 to kill B(JROW,JROW-1)

            TEMP = B( JROW, JROW );
            slartg(TEMP, B( JROW, JROW-1 ), C, S, B( JROW, JROW ) );
            B[JROW][JROW-1] = ZERO;
            srot(IHI, A( 1, JROW ), 1, A( 1, JROW-1 ), 1, C, S );
            srot(JROW-1, B( 1, JROW ), 1, B( 1, JROW-1 ), 1, C, S )             IF( ILZ ) CALL SROT( N, Z( 1, JROW ), 1, Z( 1, JROW-1 ), 1, C, S );
         } // 30
      } // 40

      }
