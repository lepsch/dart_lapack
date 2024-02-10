import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dgttrs(TRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB, Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             TRANS;
      int                INFO, LDB, N, NRHS;
      int                IPIV( * );
      double             B( LDB, * ), D( * ), DL( * ), DU( * ), DU2( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               NOTRAN;
      int                ITRANS, J, JB, NB;
      // ..
      // .. External Functions ..
      //- int                ILAENV;
      // EXTERNAL ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGTTS2, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

      INFO = 0;
      NOTRAN = ( TRANS == 'N' || TRANS == 'n' );
      if ( !NOTRAN && !( TRANS == 'T' || TRANS == 't' ) && !( TRANS == 'C' || TRANS == 'c' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( NRHS < 0 ) {
         INFO = -3;
      } else if ( LDB < max( N, 1 ) ) {
         INFO = -10;
      }
      if ( INFO != 0 ) {
         xerbla('DGTTRS', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0 || NRHS == 0) return;

      // Decode TRANS

      if ( NOTRAN ) {
         ITRANS = 0;
      } else {
         ITRANS = 1;
      }

      // Determine the number of right-hand sides to solve at a time.

      if ( NRHS == 1 ) {
         NB = 1;
      } else {
         NB = max( 1, ilaenv( 1, 'DGTTRS', TRANS, N, NRHS, -1, -1 ) );
      }

      if ( NB >= NRHS ) {
         dgtts2(ITRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB );
      } else {
         for (J = 1; NB < 0 ? J >= NRHS : J <= NRHS; J += NB) { // 10
            JB = min( NRHS-J+1, NB );
            dgtts2(ITRANS, N, JB, DL, D, DU, DU2, IPIV, B( 1, J ), LDB );
         } // 10
      }
      }
