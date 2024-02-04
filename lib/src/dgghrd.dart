      import 'dart:math';

import 'package:lapack/src/blas/drot.dart';
import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/dhgeqz.dart';
import 'package:lapack/src/dlartg.dart';
import 'package:lapack/src/dlaset.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

void dgghrd(final String COMPQ, final String COMPZ,
      final int N, final int ILO, final int IHI,
      final Matrix<double> A, final int LDA,
      final Matrix<double> B, final int LDB,
      final Matrix<double> Q, final int LDQ,
      final Matrix<double> Z, final int LDZ, final Box<int> INFO, ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      const              ONE = 1.0, ZERO = 0.0 ;
      bool               ILQ=false, ILZ=false;
      int                ICOMPQ, ICOMPZ, JCOL, JROW;
      double             C, S, TEMP;

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

      INFO.value = 0;
      if ( ICOMPQ <= 0 ) {
         INFO.value = -1;
      } else if ( ICOMPZ <= 0 ) {
         INFO.value = -2;
      } else if ( N < 0 ) {
         INFO.value = -3;
      } else if ( ILO < 1 ) {
         INFO.value = -4;
      } else if ( IHI > N || IHI < ILO-1 ) {
         INFO.value = -5;
      } else if ( LDA < max( 1, N ) ) {
         INFO.value = -7;
      } else if ( LDB < max( 1, N ) ) {
         INFO.value = -9;
      } else if ( ( ILQ && LDQ < N ) || LDQ < 1 ) {
         INFO.value = -11;
      } else if ( ( ILZ && LDZ < N ) || LDZ < 1 ) {
         INFO.value = -13;
      }
      if ( INFO.value != 0 ) {
         xerbla('DGGHRD', -INFO.value );
         return;
      }

      // Initialize Q and Z if desired.

      if (ICOMPQ == 3) dlaset( 'Full', N, N, ZERO, ONE, Q, LDQ );
      if( ICOMPZ == 3 ) dlaset( 'Full', N, N, ZERO, ONE, Z, LDZ );

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

            TEMP = A[JROW-1][ JCOL] ;
            dlartg(TEMP, A( JROW, JCOL ), C, S, A( JROW-1, JCOL ) );
            A[JROW][JCOL] = ZERO;
            drot(N-JCOL, A( JROW-1, JCOL+1 ), LDA, A( JROW, JCOL+1 ), LDA, C, S );
            drot(N+2-JROW, B( JROW-1, JROW-1 ), LDB, B( JROW, JROW-1 ), LDB, C, S );
            if( ILQ ) drot( N, Q( 1, JROW-1 ), 1, Q( 1, JROW ), 1, C, S );

            // Step 2: rotate columns JROW, JROW-1 to kill B(JROW,JROW-1)

            TEMP = B( JROW, JROW );
            dlartg(TEMP, B( JROW, JROW-1 ), C, S, B( JROW, JROW ) );
            B[JROW, JROW-1] = ZERO;
            drot(IHI, A( 1, JROW ), 1, A( 1, JROW-1 ), 1, C, S );
            drot(JROW-1, B( 1, JROW ), 1, B( 1, JROW-1 ), 1, C, S );
            if( ILZ ) drot( N, Z( 1, JROW ), 1, Z( 1, JROW-1 ), 1, C, S );
         } // 30
      } // 40
      }
