import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void drscl(N, SA, SX, INCX ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INCX, N;
      double             SA;
      double             SX( * );
      // ..

      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      bool               DONE;
      double             BIGNUM, CDEN, CDEN1, CNUM, CNUM1, MUL, SMLNUM;
      // ..
      // .. External Functions ..
      //- double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSCAL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS

      // Quick return if possible

      if (N <= 0) return;

      // Get machine parameters

      SMLNUM = dlamch( 'S' );
      BIGNUM = ONE / SMLNUM;

      // Initialize the denominator to SA and the numerator to 1.

      CDEN = SA;
      CNUM = ONE;

      } // 10
      CDEN1 = CDEN*SMLNUM;
      CNUM1 = CNUM / BIGNUM;
      if ( ( CDEN1 ).abs() > ( CNUM ).abs() && CNUM != ZERO ) {

         // Pre-multiply X by SMLNUM if CDEN is large compared to CNUM.

         MUL = SMLNUM;
         DONE = false;
         CDEN = CDEN1;
      } else if ( ( CNUM1 ).abs() > ( CDEN ).abs() ) {

         // Pre-multiply X by BIGNUM if CDEN is small compared to CNUM.

         MUL = BIGNUM;
         DONE = false;
         CNUM = CNUM1;
      } else {

         // Multiply X by CNUM / CDEN and return.

         MUL = CNUM / CDEN;
         DONE = true;
      }

      // Scale the vector X by MUL

      dscal(N, MUL, SX, INCX );

      if ( !DONE) GO TO 10;

      return;
      }
