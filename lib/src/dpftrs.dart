import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dpftrs(TRANSR, UPLO, N, NRHS, A, B, LDB, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             TRANSR, UPLO;
      int                INFO, LDB, N, NRHS;
      double             A( 0: * ), B( LDB, * );
      // ..

      double             ONE;
      const              ONE = 1.0 ;
      bool               LOWER, NORMALTRANSR;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, DTFSM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      // Test the input parameters.

      INFO = 0;
      NORMALTRANSR = lsame( TRANSR, 'N' );
      LOWER = lsame( UPLO, 'L' );
      if ( !NORMALTRANSR && !lsame( TRANSR, 'T' ) ) {
         INFO = -1;
      } else if ( !LOWER && !lsame( UPLO, 'U' ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( NRHS < 0 ) {
         INFO = -4;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -7;
      }
      if ( INFO != 0 ) {
         xerbla('DPFTRS', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0 || NRHS == 0) return;

      // start execution: there are two triangular solves

      if ( LOWER ) {
         dtfsm(TRANSR, 'L', UPLO, 'N', 'N', N, NRHS, ONE, A, B, LDB );
         dtfsm(TRANSR, 'L', UPLO, 'T', 'N', N, NRHS, ONE, A, B, LDB );
      } else {
         dtfsm(TRANSR, 'L', UPLO, 'T', 'N', N, NRHS, ONE, A, B, LDB );
         dtfsm(TRANSR, 'L', UPLO, 'N', 'N', N, NRHS, ONE, A, B, LDB );
      }

      return;
      }
