import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dlaqz1(A, LDA, B, LDB, SR1, SR2, SI, BETA1, BETA2, V ) {
      // IMPLICIT NONE

      // Arguments
      int    , INTENT( IN ) :: LDA, LDB;
      double          , INTENT( IN ) :: A( LDA, * ), B( LDB, * ), SR1, SR2, SI, BETA1, BETA2;
      double          , INTENT( OUT ) :: V( * );

      // Parameters
      double           :: ZERO, ONE, HALF;
      const    ZERO = 0.0, ONE = 1.0, HALF = 0.5 ;

      // Local scalars
      double           :: W( 2 ), SAFMIN, SAFMAX, SCALE1, SCALE2;

      // External Functions
      //- double          , EXTERNAL :: DLAMCH;
      //- bool   , EXTERNAL :: DISNAN;

      SAFMIN = dlamch( 'SAFE MINIMUM' );
      SAFMAX = ONE/SAFMIN;

      // Calculate first shifted vector

      W[1] = BETA1*A( 1, 1 )-SR1*B( 1, 1 );
      W[2] = BETA1*A( 2, 1 )-SR1*B( 2, 1 );
      SCALE1 = sqrt( ( W( 1 ) ).abs() ) * sqrt( ( W( 2 ) ).abs() );
      if ( SCALE1 >= SAFMIN && SCALE1 <= SAFMAX ) {
         W[1] = W( 1 )/SCALE1;
         W[2] = W( 2 )/SCALE1;
      }

      // Solve linear system

      W[2] = W( 2 )/B( 2, 2 );
      W[1] = ( W( 1 )-B( 1, 2 )*W( 2 ) )/B( 1, 1 );
      SCALE2 = sqrt( ( W( 1 ) ).abs() ) * sqrt( ( W( 2 ) ).abs() );
      if ( SCALE2 >= SAFMIN && SCALE2 <= SAFMAX ) {
         W[1] = W( 1 )/SCALE2;
         W[2] = W( 2 )/SCALE2;
      }

      // Apply second shift

      V[1] = BETA2*( A( 1, 1 )*W( 1 )+A( 1, 2 )*W( 2 ) )-SR2*( B( 1, 1 )*W( 1 )+B( 1, 2 )*W( 2 ) )       V( 2 ) = BETA2*( A( 2, 1 )*W( 1 )+A( 2, 2 )*W( 2 ) )-SR2*( B( 2, 1 )*W( 1 )+B( 2, 2 )*W( 2 ) )       V( 3 ) = BETA2*( A( 3, 1 )*W( 1 )+A( 3, 2 )*W( 2 ) )-SR2*( B( 3, 1 )*W( 1 )+B( 3, 2 )*W( 2 ) );

      // Account for imaginary part

      V[1] = V( 1 )+SI*SI*B( 1, 1 )/SCALE1/SCALE2;

      // Check for overflow

      if ( ( V( 1 ) ).abs() > SAFMAX || ( V( 2 ) ).abs() > SAFMAX || ( V( 3 ) ).abs() > SAFMAX || disnan( V( 1 ) ) || disnan( V( 2 ) ) || disnan( V( 3 ) ) ) {
         V[1] = ZERO;
         V[2] = ZERO;
         V[3] = ZERO;
      }
      END SUBROUTINE;
