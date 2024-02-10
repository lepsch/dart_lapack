import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dtprfs(UPLO, TRANS, DIAG, N, NRHS, AP, B, LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             DIAG, TRANS, UPLO;
      int                INFO, LDB, LDX, N, NRHS;
      int                IWORK( * );
      double             AP( * ), B( LDB, * ), BERR( * ), FERR( * ), WORK( * ), X( LDX, * );
      // ..

      double             ZERO;
      const              ZERO = 0.0 ;
      double             ONE;
      const              ONE = 1.0 ;
      bool               NOTRAN, NOUNIT, UPPER;
      String             TRANST;
      int                I, J, K, KASE, KC, NZ;
      double             EPS, LSTRES, S, SAFE1, SAFE2, SAFMIN, XK;
      int                ISAVE( 3 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DCOPY, DLACN2, DTPMV, DTPSV, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DLAMCH;
      // EXTERNAL lsame, DLAMCH

      // Test the input parameters.

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      NOTRAN = lsame( TRANS, 'N' );
      NOUNIT = lsame( DIAG, 'N' );

      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( !NOTRAN && !lsame( TRANS, 'T' ) && !lsame( TRANS, 'C' ) ) {
         INFO = -2;
      } else if ( !NOUNIT && !lsame( DIAG, 'U' ) ) {
         INFO = -3;
      } else if ( N < 0 ) {
         INFO = -4;
      } else if ( NRHS < 0 ) {
         INFO = -5;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -8;
      } else if ( LDX < max( 1, N ) ) {
         INFO = -10;
      }
      if ( INFO != 0 ) {
         xerbla('DTPRFS', -INFO );
         return;
      }

      // Quick return if possible

      if ( N == 0 || NRHS == 0 ) {
         for (J = 1; J <= NRHS; J++) { // 10
            FERR[J] = ZERO;
            BERR[J] = ZERO;
         } // 10
         return;
      }

      if ( NOTRAN ) {
         TRANST = 'T';
      } else {
         TRANST = 'N';
      }

      // NZ = maximum number of nonzero elements in each row of A, plus 1

      NZ = N + 1;
      EPS = dlamch( 'Epsilon' );
      SAFMIN = dlamch( 'Safe minimum' );
      SAFE1 = NZ*SAFMIN;
      SAFE2 = SAFE1 / EPS;

      // Do for each right hand side

      for (J = 1; J <= NRHS; J++) { // 250

         // Compute residual R = B - op(A) * X,
         // where op(A) = A or A**T, depending on TRANS.

         dcopy(N, X( 1, J ), 1, WORK( N+1 ), 1 );
         dtpmv(UPLO, TRANS, DIAG, N, AP, WORK( N+1 ), 1 );
         daxpy(N, -ONE, B( 1, J ), 1, WORK( N+1 ), 1 );

         // Compute componentwise relative backward error from formula

         // max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) )

         // where abs(Z) is the componentwise absolute value of the matrix
         // or vector Z.  If the i-th component of the denominator is less
         // than SAFE2, then SAFE1 is added to the i-th components of the
         // numerator and denominator before dividing.

         for (I = 1; I <= N; I++) { // 20
            WORK[I] = ( B( I, J ) ).abs();
         } // 20

         if ( NOTRAN ) {

            // Compute abs(A)*abs(X) + abs(B).

            if ( UPPER ) {
               KC = 1;
               if ( NOUNIT ) {
                  for (K = 1; K <= N; K++) { // 40
                     XK = ( X( K, J ) ).abs();
                     for (I = 1; I <= K; I++) { // 30
                        WORK[I] = WORK( I ) + ( AP( KC+I-1 ) ).abs()*XK;
                     } // 30
                     KC = KC + K;
                  } // 40
               } else {
                  for (K = 1; K <= N; K++) { // 60
                     XK = ( X( K, J ) ).abs();
                     for (I = 1; I <= K - 1; I++) { // 50
                        WORK[I] = WORK( I ) + ( AP( KC+I-1 ) ).abs()*XK;
                     } // 50
                     WORK[K] = WORK( K ) + XK;
                     KC = KC + K;
                  } // 60
               }
            } else {
               KC = 1;
               if ( NOUNIT ) {
                  for (K = 1; K <= N; K++) { // 80
                     XK = ( X( K, J ) ).abs();
                     for (I = K; I <= N; I++) { // 70
                        WORK[I] = WORK( I ) + ( AP( KC+I-K ) ).abs()*XK;
                     } // 70
                     KC = KC + N - K + 1;
                  } // 80
               } else {
                  for (K = 1; K <= N; K++) { // 100
                     XK = ( X( K, J ) ).abs();
                     for (I = K + 1; I <= N; I++) { // 90
                        WORK[I] = WORK( I ) + ( AP( KC+I-K ) ).abs()*XK;
                     } // 90
                     WORK[K] = WORK( K ) + XK;
                     KC = KC + N - K + 1;
                  } // 100
               }
            }
         } else {

            // Compute abs(A**T)*abs(X) + abs(B).

            if ( UPPER ) {
               KC = 1;
               if ( NOUNIT ) {
                  for (K = 1; K <= N; K++) { // 120
                     S = ZERO;
                     for (I = 1; I <= K; I++) { // 110
                        S = S + ( AP( KC+I-1 ) ).abs()*( X( I, J ) ).abs();
                     } // 110
                     WORK[K] = WORK( K ) + S;
                     KC = KC + K;
                  } // 120
               } else {
                  for (K = 1; K <= N; K++) { // 140
                     S = ( X( K, J ) ).abs();
                     for (I = 1; I <= K - 1; I++) { // 130
                        S = S + ( AP( KC+I-1 ) ).abs()*( X( I, J ) ).abs();
                     } // 130
                     WORK[K] = WORK( K ) + S;
                     KC = KC + K;
                  } // 140
               }
            } else {
               KC = 1;
               if ( NOUNIT ) {
                  for (K = 1; K <= N; K++) { // 160
                     S = ZERO;
                     for (I = K; I <= N; I++) { // 150
                        S = S + ( AP( KC+I-K ) ).abs()*( X( I, J ) ).abs();
                     } // 150
                     WORK[K] = WORK( K ) + S;
                     KC = KC + N - K + 1;
                  } // 160
               } else {
                  for (K = 1; K <= N; K++) { // 180
                     S = ( X( K, J ) ).abs();
                     for (I = K + 1; I <= N; I++) { // 170
                        S = S + ( AP( KC+I-K ) ).abs()*( X( I, J ) ).abs();
                     } // 170
                     WORK[K] = WORK( K ) + S;
                     KC = KC + N - K + 1;
                  } // 180
               }
            }
         }
         S = ZERO;
         for (I = 1; I <= N; I++) { // 190
            if ( WORK( I ) > SAFE2 ) {
               S = max( S, ( WORK( N+I ) ).abs() / WORK( I ) );
            } else {
               S = max( S, ( ( WORK( N+I ) ).abs()+SAFE1 ) / ( WORK( I )+SAFE1 ) );
            }
         } // 190
         BERR[J] = S;

         // Bound error from formula

         // norm(X - XTRUE) / norm(X) <= FERR =
         // norm( abs(inv(op(A)))*
         //    ( abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) ))) / norm(X)

         // where
         //   norm(Z) is the magnitude of the largest component of Z
         //   inv(op(A)) is the inverse of op(A)
         //   abs(Z) is the componentwise absolute value of the matrix or
         //      vector Z
         //   NZ is the maximum number of nonzeros in any row of A, plus 1
         //   EPS is machine epsilon

         // The i-th component of abs(R)+NZ*EPS*(abs(op(A))*abs(X)+abs(B))
         // is incremented by SAFE1 if the i-th component of
         // abs(op(A))*abs(X) + abs(B) is less than SAFE2.

         // Use DLACN2 to estimate the infinity-norm of the matrix
         //    inv(op(A)) * diag(W),
         // where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) )))

         for (I = 1; I <= N; I++) { // 200
            if ( WORK( I ) > SAFE2 ) {
               WORK[I] = ( WORK( N+I ) ).abs() + NZ*EPS*WORK( I );
            } else {
               WORK[I] = ( WORK( N+I ) ).abs() + NZ*EPS*WORK( I ) + SAFE1;
            }
         } // 200

         KASE = 0;
         } // 210
         dlacn2(N, WORK( 2*N+1 ), WORK( N+1 ), IWORK, FERR( J ), KASE, ISAVE );
         if ( KASE != 0 ) {
            if ( KASE == 1 ) {

               // Multiply by diag(W)*inv(op(A)**T).

               dtpsv(UPLO, TRANST, DIAG, N, AP, WORK( N+1 ), 1 );
               for (I = 1; I <= N; I++) { // 220
                  WORK[N+I] = WORK( I )*WORK( N+I );
               } // 220
            } else {

               // Multiply by inv(op(A))*diag(W).

               for (I = 1; I <= N; I++) { // 230
                  WORK[N+I] = WORK( I )*WORK( N+I );
               } // 230
               dtpsv(UPLO, TRANS, DIAG, N, AP, WORK( N+1 ), 1 );
            }
            GO TO 210;
         }

         // Normalize error.

         LSTRES = ZERO;
         for (I = 1; I <= N; I++) { // 240
            LSTRES = max( LSTRES, ( X( I, J ) ).abs() );
         } // 240
         if (LSTRES != ZERO) FERR( J ) = FERR( J ) / LSTRES;

      } // 250

      }
