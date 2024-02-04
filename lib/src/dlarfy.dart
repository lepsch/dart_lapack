import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dlarfy(UPLO, N, V, INCV, TAU, C, LDC, WORK ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INCV, LDC, N;
      double             TAU;
      // ..
      // .. Array Arguments ..
      double             C( LDC, * ), V( * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE, ZERO, HALF;
      const              ONE = 1.0, ZERO = 0.0, HALF = 0.5 ;
      // ..
      // .. Local Scalars ..
      double             ALPHA;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DSYMV, DSYR2
      // ..
      // .. External Functions ..
      //- double             DDOT;
      // EXTERNAL DDOT
      // ..
      // .. Executable Statements ..

      if (TAU == ZERO) return;

      // Form  w:= C * v

      dsymv(UPLO, N, ONE, C, LDC, V, INCV, ZERO, WORK, 1 );

      ALPHA = -HALF*TAU*ddot( N, WORK, 1, V, INCV );
      daxpy(N, ALPHA, V, INCV, WORK, 1 );

      // C := C - v * w' - w * v'

      dsyr2(UPLO, N, -TAU, V, INCV, WORK, 1, C, LDC );

      return;
      }
