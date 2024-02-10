import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      RECURSIVE SUBROUTINE DPOTRF2( UPLO, N, A, LDA, INFO );

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, LDA, N;
      double             A( LDA, * );
      // ..

      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      bool               UPPER;
      int                N1, N2, IINFO;
      // ..
      // .. External Functions ..
      //- bool               lsame, DISNAN;
      // EXTERNAL lsame, DISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSYRK, DTRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT

      // Test the input parameters

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -4;
      }
      if ( INFO != 0 ) {
         xerbla('DPOTRF2', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // N=1 case

      if ( N == 1 ) {

         // Test for non-positive-definiteness

         if ( A( 1, 1 ) <= ZERO || disnan( A( 1, 1 ) ) ) {
            INFO = 1;
            return;
         }

         // Factor

         A[1][1] = sqrt( A( 1, 1 ) );

      // Use recursive code

      } else {
         N1 = N/2;
         N2 = N-N1;

         // Factor A11

         dpotrf2(UPLO, N1, A( 1, 1 ), LDA, IINFO );
         if ( IINFO != 0 ) {
            INFO = IINFO;
            return;
         }

         // Compute the Cholesky factorization A = U**T*U

         if ( UPPER ) {

            // Update and scale A12

            dtrsm('L', 'U', 'T', 'N', N1, N2, ONE, A( 1, 1 ), LDA, A( 1, N1+1 ), LDA );

            // Update and factor A22

            dsyrk(UPLO, 'T', N2, N1, -ONE, A( 1, N1+1 ), LDA, ONE, A( N1+1, N1+1 ), LDA );
            dpotrf2(UPLO, N2, A( N1+1, N1+1 ), LDA, IINFO );
            if ( IINFO != 0 ) {
               INFO = IINFO + N1;
               return;
            }

         // Compute the Cholesky factorization A = L*L**T

         } else {

            // Update and scale A21

            dtrsm('R', 'L', 'T', 'N', N2, N1, ONE, A( 1, 1 ), LDA, A( N1+1, 1 ), LDA );

            // Update and factor A22

            dsyrk(UPLO, 'N', N2, N1, -ONE, A( N1+1, 1 ), LDA, ONE, A( N1+1, N1+1 ), LDA );
            dpotrf2(UPLO, N2, A( N1+1, N1+1 ), LDA, IINFO );
            if ( IINFO != 0 ) {
               INFO = IINFO + N1;
               return;
            }
         }
      }
      }
