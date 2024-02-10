import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dporfs(UPLO, N, NRHS, final Matrix<double> A, final int LDA, final Matrix<double> AF, final int LDAF, final Matrix<double> B, final int LDB, final Matrix<double> X, final int LDX, FERR, BERR, WORK, IWORK, Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, LDA, LDAF, LDB, LDX, N, NRHS;
      int                IWORK( * );
      double             A( LDA, * ), AF( LDAF, * ), B( LDB, * ), BERR( * ), FERR( * ), WORK( * ), X( LDX, * );
      // ..

      int                ITMAX;
      const              ITMAX = 5 ;
      double             ZERO;
      const              ZERO = 0.0 ;
      double             ONE;
      const              ONE = 1.0 ;
      double             TWO;
      const              TWO = 2.0 ;
      double             THREE;
      const              THREE = 3.0 ;
      bool               UPPER;
      int                COUNT, I, J, K, KASE, NZ;
      double             EPS, LSTRES, S, SAFE1, SAFE2, SAFMIN, XK;
      int                ISAVE( 3 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DCOPY, DLACN2, DPOTRS, DSYMV, XERBLA
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
      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( NRHS < 0 ) {
         INFO = -3;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -5;
      } else if ( LDAF < max( 1, N ) ) {
         INFO = -7;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -9;
      } else if ( LDX < max( 1, N ) ) {
         INFO = -11;
      }
      if ( INFO != 0 ) {
         xerbla('DPORFS', -INFO );
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

      // NZ = maximum number of nonzero elements in each row of A, plus 1

      NZ = N + 1;
      EPS = dlamch( 'Epsilon' );
      SAFMIN = dlamch( 'Safe minimum' );
      SAFE1 = NZ*SAFMIN;
      SAFE2 = SAFE1 / EPS;

      // Do for each right hand side

      for (J = 1; J <= NRHS; J++) { // 140

         COUNT = 1;
         LSTRES = THREE;
         } // 20

         // Loop until stopping criterion is satisfied.

         // Compute residual R = B - A * X

         dcopy(N, B( 1, J ), 1, WORK( N+1 ), 1 );
         dsymv(UPLO, N, -ONE, A, LDA, X( 1, J ), 1, ONE, WORK( N+1 ), 1 );

         // Compute componentwise relative backward error from formula

         // max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) )

         // where abs(Z) is the componentwise absolute value of the matrix
         // or vector Z.  If the i-th component of the denominator is less
         // than SAFE2, then SAFE1 is added to the i-th components of the
         // numerator and denominator before dividing.

         for (I = 1; I <= N; I++) { // 30
            WORK[I] = ( B( I, J ) ).abs();
         } // 30

         // Compute abs(A)*abs(X) + abs(B).

         if ( UPPER ) {
            for (K = 1; K <= N; K++) { // 50
               S = ZERO;
               XK = ( X( K, J ) ).abs();
               for (I = 1; I <= K - 1; I++) { // 40
                  WORK[I] = WORK( I ) + ( A( I, K ) ).abs()*XK;
                  S = S + ( A( I, K ) ).abs()*( X( I, J ) ).abs();
               } // 40
               WORK[K] = WORK( K ) + ( A( K, K ) ).abs()*XK + S;
            } // 50
         } else {
            for (K = 1; K <= N; K++) { // 70
               S = ZERO;
               XK = ( X( K, J ) ).abs();
               WORK[K] = WORK( K ) + ( A( K, K ) ).abs()*XK;
               for (I = K + 1; I <= N; I++) { // 60
                  WORK[I] = WORK( I ) + ( A( I, K ) ).abs()*XK;
                  S = S + ( A( I, K ) ).abs()*( X( I, J ) ).abs();
               } // 60
               WORK[K] = WORK( K ) + S;
            } // 70
         }
         S = ZERO;
         for (I = 1; I <= N; I++) { // 80
            if ( WORK( I ) > SAFE2 ) {
               S = max( S, ( WORK( N+I ) ).abs() / WORK( I ) );
            } else {
               S = max( S, ( ( WORK( N+I ) ).abs()+SAFE1 ) / ( WORK( I )+SAFE1 ) );
            }
         } // 80
         BERR[J] = S;

         // Test stopping criterion. Continue iterating if
         //    1) The residual BERR(J) is larger than machine epsilon, and
         //    2) BERR(J) decreased by at least a factor of 2 during the
         //       last iteration, and
         //    3) At most ITMAX iterations tried.

         if ( BERR( J ) > EPS && TWO*BERR( J ) <= LSTRES && COUNT <= ITMAX ) {

            // Update solution and try again.

            dpotrs(UPLO, N, 1, AF, LDAF, WORK( N+1 ), N, INFO );
            daxpy(N, ONE, WORK( N+1 ), 1, X( 1, J ), 1 );
            LSTRES = BERR( J );
            COUNT = COUNT + 1;
            GO TO 20;
         }

         // Bound error from formula

         // norm(X - XTRUE) / norm(X) <= FERR =
         // norm( abs(inv(A))*
         //    ( abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) ))) / norm(X)

         // where
         //   norm(Z) is the magnitude of the largest component of Z
         //   inv(A) is the inverse of A
         //   abs(Z) is the componentwise absolute value of the matrix or
         //      vector Z
         //   NZ is the maximum number of nonzeros in any row of A, plus 1
         //   EPS is machine epsilon

         // The i-th component of abs(R)+NZ*EPS*(abs(A)*abs(X)+abs(B))
         // is incremented by SAFE1 if the i-th component of
         // abs(A)*abs(X) + abs(B) is less than SAFE2.

         // Use DLACN2 to estimate the infinity-norm of the matrix
         //    inv(A) * diag(W),
         // where W = abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) )))

         for (I = 1; I <= N; I++) { // 90
            if ( WORK( I ) > SAFE2 ) {
               WORK[I] = ( WORK( N+I ) ).abs() + NZ*EPS*WORK( I );
            } else {
               WORK[I] = ( WORK( N+I ) ).abs() + NZ*EPS*WORK( I ) + SAFE1;
            }
         } // 90

         KASE = 0;
         } // 100
         dlacn2(N, WORK( 2*N+1 ), WORK( N+1 ), IWORK, FERR( J ), KASE, ISAVE );
         if ( KASE != 0 ) {
            if ( KASE == 1 ) {

               // Multiply by diag(W)*inv(A**T).

               dpotrs(UPLO, N, 1, AF, LDAF, WORK( N+1 ), N, INFO );
               for (I = 1; I <= N; I++) { // 110
                  WORK[N+I] = WORK( I )*WORK( N+I );
               } // 110
            } else if ( KASE == 2 ) {

               // Multiply by inv(A)*diag(W).

               for (I = 1; I <= N; I++) { // 120
                  WORK[N+I] = WORK( I )*WORK( N+I );
               } // 120
               dpotrs(UPLO, N, 1, AF, LDAF, WORK( N+1 ), N, INFO );
            }
            GO TO 100;
         }

         // Normalize error.

         LSTRES = ZERO;
         for (I = 1; I <= N; I++) { // 130
            LSTRES = max( LSTRES, ( X( I, J ) ).abs() );
         } // 130
         if (LSTRES != ZERO) FERR( J ) = FERR( J ) / LSTRES;

      } // 140

      }
