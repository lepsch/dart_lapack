import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dlarz(SIDE, M, N, L, V, INCV, TAU, final Matrix<double> C, final int LDC, WORK ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             SIDE;
      int                INCV, L, LDC, M, N;
      double             TAU;
      double             C( LDC, * ), V( * ), WORK( * );
      // ..

      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DCOPY, DGEMV, DGER
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame

      if ( lsame( SIDE, 'L' ) ) {

         // Form  H * C

         if ( TAU != ZERO ) {

            // w( 1:n ) = C( 1, 1:n )

            dcopy(N, C, LDC, WORK, 1 );

            // w( 1:n ) = w( 1:n ) + C( m-l+1:m, 1:n )**T * v( 1:l )

            dgemv('Transpose', L, N, ONE, C( M-L+1, 1 ), LDC, V, INCV, ONE, WORK, 1 );

            // C( 1, 1:n ) = C( 1, 1:n ) - tau * w( 1:n )

            daxpy(N, -TAU, WORK, 1, C, LDC );

            // C( m-l+1:m, 1:n ) = C( m-l+1:m, 1:n ) - ...
            //                     tau * v( 1:l ) * w( 1:n )**T

            dger(L, N, -TAU, V, INCV, WORK, 1, C( M-L+1, 1 ), LDC );
         }

      } else {

         // Form  C * H

         if ( TAU != ZERO ) {

            // w( 1:m ) = C( 1:m, 1 )

            dcopy(M, C, 1, WORK, 1 );

            // w( 1:m ) = w( 1:m ) + C( 1:m, n-l+1:n, 1:n ) * v( 1:l )

            dgemv('No transpose', M, L, ONE, C( 1, N-L+1 ), LDC, V, INCV, ONE, WORK, 1 );

            // C( 1:m, 1 ) = C( 1:m, 1 ) - tau * w( 1:m )

            daxpy(M, -TAU, WORK, 1, C, 1 );

            // C( 1:m, n-l+1:n ) = C( 1:m, n-l+1:n ) - ...
            //                     tau * w( 1:m ) * v( 1:l )**T

            dger(M, L, -TAU, WORK, 1, V, INCV, C( 1, N-L+1 ), LDC );

         }

      }

      }
