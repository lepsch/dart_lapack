import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dptrfs(N, NRHS, D, E, DF, EF, final Matrix<double> B, final int LDB, final Matrix<double> X, final int LDX, FERR, BERR, WORK, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDB, LDX, N, NRHS;
      double             B( LDB, * ), BERR( * ), D( * ), DF( * ), E( * ), EF( * ), FERR( * ), WORK( * ), X( LDX, * );
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
      int                COUNT, I, IX, J, NZ;
      double             BI, CX, DX, EPS, EX, LSTRES, S, SAFE1, SAFE2, SAFMIN;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DPTTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. External Functions ..
      //- int                idamax;
      //- double             DLAMCH;
      // EXTERNAL idamax, DLAMCH

      // Test the input parameters.

      INFO = 0;
      if ( N < 0 ) {
         INFO = -1;
      } else if ( NRHS < 0 ) {
         INFO = -2;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -8;
      } else if ( LDX < max( 1, N ) ) {
         INFO = -10;
      }
      if ( INFO != 0 ) {
         xerbla('DPTRFS', -INFO );
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

      NZ = 4;
      EPS = dlamch( 'Epsilon' );
      SAFMIN = dlamch( 'Safe minimum' );
      SAFE1 = NZ*SAFMIN;
      SAFE2 = SAFE1 / EPS;

      // Do for each right hand side

      for (J = 1; J <= NRHS; J++) { // 90

         COUNT = 1;
         LSTRES = THREE;
         } // 20

         // Loop until stopping criterion is satisfied.

         // Compute residual R = B - A * X.  Also compute
         // abs(A)*abs(x) + abs(b) for use in the backward error bound.

         if ( N == 1 ) {
            BI = B( 1, J );
            DX = D( 1 )*X( 1, J );
            WORK[N+1] = BI - DX;
            WORK[1] = ( BI ).abs() + ( DX ).abs();
         } else {
            BI = B( 1, J );
            DX = D( 1 )*X( 1, J );
            EX = E( 1 )*X( 2, J );
            WORK[N+1] = BI - DX - EX;
            WORK[1] = ( BI ).abs() + ( DX ).abs() + ( EX ).abs();
            for (I = 2; I <= N - 1; I++) { // 30
               BI = B( I, J );
               CX = E( I-1 )*X( I-1, J );
               DX = D( I )*X( I, J );
               EX = E( I )*X( I+1, J );
               WORK[N+I] = BI - CX - DX - EX;
               WORK[I] = ( BI ).abs() + ( CX ).abs() + ( DX ).abs() + ( EX ).abs();
            } // 30
            BI = B( N, J );
            CX = E( N-1 )*X( N-1, J );
            DX = D( N )*X( N, J );
            WORK[N+N] = BI - CX - DX;
            WORK[N] = ( BI ).abs() + ( CX ).abs() + ( DX ).abs();
         }

         // Compute componentwise relative backward error from formula

         // max(i) ( abs(R(i)) / ( abs(A)*abs(X) + abs(B) )(i) )

         // where abs(Z) is the componentwise absolute value of the matrix
         // or vector Z.  If the i-th component of the denominator is less
         // than SAFE2, then SAFE1 is added to the i-th components of the
         // numerator and denominator before dividing.

         S = ZERO;
         for (I = 1; I <= N; I++) { // 40
            if ( WORK( I ) > SAFE2 ) {
               S = max( S, ( WORK( N+I ) ).abs() / WORK( I ) );
            } else {
               S = max( S, ( ( WORK( N+I ) ).abs()+SAFE1 ) / ( WORK( I )+SAFE1 ) );
            }
         } // 40
         BERR[J] = S;

         // Test stopping criterion. Continue iterating if
         //    1) The residual BERR(J) is larger than machine epsilon, and
         //    2) BERR(J) decreased by at least a factor of 2 during the
         //       last iteration, and
         //    3) At most ITMAX iterations tried.

         if ( BERR( J ) > EPS && TWO*BERR( J ) <= LSTRES && COUNT <= ITMAX ) {

            // Update solution and try again.

            dpttrs(N, 1, DF, EF, WORK( N+1 ), N, INFO );
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

         for (I = 1; I <= N; I++) { // 50
            if ( WORK( I ) > SAFE2 ) {
               WORK[I] = ( WORK( N+I ) ).abs() + NZ*EPS*WORK( I );
            } else {
               WORK[I] = ( WORK( N+I ) ).abs() + NZ*EPS*WORK( I ) + SAFE1;
            }
         } // 50
         IX = idamax( N, WORK, 1 );
         FERR[J] = WORK( IX );

         // Estimate the norm of inv(A).

         // Solve M(A) * x = e, where M(A) = (m(i,j)) is given by

            // m(i,j) =  abs(A(i,j)), i = j,
            // m(i,j) = -abs(A(i,j)), i != j,

         // and e = [ 1, 1, ..., 1 ]**T.  Note M(A) = M(L)*D*M(L)**T.

         // Solve M(L) * x = e.

         WORK[1] = ONE;
         for (I = 2; I <= N; I++) { // 60
            WORK[I] = ONE + WORK( I-1 )*( EF( I-1 ) ).abs();
         } // 60

         // Solve D * M(L)**T * x = b.

         WORK[N] = WORK( N ) / DF( N );
         for (I = N - 1; I >= 1; I--) { // 70
            WORK[I] = WORK( I ) / DF( I ) + WORK( I+1 )*( EF( I ) ).abs();
         } // 70

         // Compute norm(inv(A)) = max(x(i)), 1<=i<=n.

         IX = idamax( N, WORK, 1 );
         FERR[J] = FERR( J )*( WORK( IX ) ).abs();

         // Normalize error.

         LSTRES = ZERO;
         for (I = 1; I <= N; I++) { // 80
            LSTRES = max( LSTRES, ( X( I, J ) ).abs() );
         } // 80
         if (LSTRES != ZERO) FERR( J ) = FERR( J ) / LSTRES;

      } // 90

      }
