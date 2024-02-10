import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dptcon(final int N, final int D, final int E, final int ANORM, final int RCOND, final Array<double> _WORK, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, N;
      double             ANORM, RCOND;
      double             D( * ), E( * ), WORK( * );
      // ..

      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                I, IX;
      double             AINVNM;
      // ..
      // .. External Functions ..
      //- int                idamax;
      // EXTERNAL idamax
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS

      // Test the input arguments.

      INFO = 0;
      if ( N < 0 ) {
         INFO = -1;
      } else if ( ANORM < ZERO ) {
         INFO = -4;
      }
      if ( INFO != 0 ) {
         xerbla('DPTCON', -INFO );
         return;
      }

      // Quick return if possible

      RCOND = ZERO;
      if ( N == 0 ) {
         RCOND = ONE;
         return;
      } else if ( ANORM == ZERO ) {
         return;
      }

      // Check that D(1:N) is positive.

      for (I = 1; I <= N; I++) { // 10
         if( D( I ) <= ZERO ) return;
      } // 10

      // Solve M(A) * x = e, where M(A) = (m(i,j)) is given by

         // m(i,j) =  abs(A(i,j)), i = j,
         // m(i,j) = -abs(A(i,j)), i != j,

      // and e = [ 1, 1, ..., 1 ]**T.  Note M(A) = M(L)*D*M(L)**T.

      // Solve M(L) * x = e.

      WORK[1] = ONE;
      for (I = 2; I <= N; I++) { // 20
         WORK[I] = ONE + WORK( I-1 )*( E( I-1 ) ).abs();
      } // 20

      // Solve D * M(L)**T * x = b.

      WORK[N] = WORK( N ) / D( N );
      for (I = N - 1; I >= 1; I--) { // 30
         WORK[I] = WORK( I ) / D( I ) + WORK( I+1 )*( E( I ) ).abs();
      } // 30

      // Compute AINVNM = max(x(i)), 1<=i<=n.

      IX = idamax( N, WORK, 1 );
      AINVNM = ( WORK( IX ) ).abs();

      // Compute the reciprocal condition number.

      if (AINVNM != ZERO) RCOND = ( ONE / AINVNM ) / ANORM;

      }
