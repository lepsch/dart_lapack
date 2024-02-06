import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dgebd2(M, N, A, LDA, D, E, TAUQ, TAUP, WORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, M, N;
      double             A( LDA, * ), D( * ), E( * ), TAUP( * ), TAUQ( * ), WORK( * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                I;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLARF, DLARFG, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

      // Test the input parameters

      INFO = 0;
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -4;
      }
      if ( INFO < 0 ) {
         xerbla('DGEBD2', -INFO );
         return;
      }

      if ( M >= N ) {

         // Reduce to upper bidiagonal form

         for (I = 1; I <= N; I++) { // 10

            // Generate elementary reflector H(i) to annihilate A(i+1:m,i)

            dlarfg(M-I+1, A( I, I ), A( min( I+1, M ), I ), 1, TAUQ( I ) );
            D[I] = A( I, I );
            A[I][I] = ONE;

            // Apply H(i) to A(i:m,i+1:n) from the left

            if (I < N) dlarf( 'Left', M-I+1, N-I, A( I, I ), 1, TAUQ( I ), A( I, I+1 ), LDA, WORK );
            A[I][I] = D( I );

            if ( I < N ) {

               // Generate elementary reflector G(i) to annihilate
               // A(i,i+2:n)

               dlarfg(N-I, A( I, I+1 ), A( I, min( I+2, N ) ), LDA, TAUP( I ) );
               E[I] = A( I, I+1 );
               A[I][I+1] = ONE;

               // Apply G(i) to A(i+1:m,i+1:n) from the right

               dlarf('Right', M-I, N-I, A( I, I+1 ), LDA, TAUP( I ), A( I+1, I+1 ), LDA, WORK );
               A[I][I+1] = E( I );
            } else {
               TAUP[I] = ZERO;
            }
         } // 10
      } else {

         // Reduce to lower bidiagonal form

         for (I = 1; I <= M; I++) { // 20

            // Generate elementary reflector G(i) to annihilate A(i,i+1:n)

            dlarfg(N-I+1, A( I, I ), A( I, min( I+1, N ) ), LDA, TAUP( I ) );
            D[I] = A( I, I );
            A[I][I] = ONE;

            // Apply G(i) to A(i+1:m,i:n) from the right

            if (I < M) dlarf( 'Right', M-I, N-I+1, A( I, I ), LDA, TAUP( I ), A( I+1, I ), LDA, WORK );
            A[I][I] = D( I );

            if ( I < M ) {

               // Generate elementary reflector H(i) to annihilate
               // A(i+2:m,i)

               dlarfg(M-I, A( I+1, I ), A( min( I+2, M ), I ), 1, TAUQ( I ) );
               E[I] = A( I+1, I );
               A[I+1][I] = ONE;

               // Apply H(i) to A(i+1:m,i+1:n) from the left

               dlarf('Left', M-I, N-I, A( I+1, I ), 1, TAUQ( I ), A( I+1, I+1 ), LDA, WORK );
               A[I+1][I] = E( I );
            } else {
               TAUQ[I] = ZERO;
            }
         } // 20
      }
      return;
      }
