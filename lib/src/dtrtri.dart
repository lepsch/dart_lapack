import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dtrtri(UPLO, DIAG, N, A, LDA, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             DIAG, UPLO;
      int                INFO, LDA, N;
      double             A( LDA, * );
      // ..

      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      bool               NOUNIT, UPPER;
      int                J, JB, NB, NN;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV;
      // EXTERNAL lsame, ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL DTRMM, DTRSM, DTRTI2, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

      // Test the input parameters.

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      NOUNIT = lsame( DIAG, 'N' );
      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( !NOUNIT && !lsame( DIAG, 'U' ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -5;
      }
      if ( INFO != 0 ) {
         xerbla('DTRTRI', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Check for singularity if non-unit.

      if ( NOUNIT ) {
         for (INFO = 1; INFO <= N; INFO++) { // 10
            if( A( INFO, INFO ) == ZERO ) return;
         } // 10
         INFO = 0;
      }

      // Determine the block size for this environment.

      NB = ilaenv( 1, 'DTRTRI', UPLO + DIAG, N, -1, -1, -1 );
      if ( NB <= 1 || NB >= N ) {

         // Use unblocked code

         dtrti2(UPLO, DIAG, N, A, LDA, INFO );
      } else {

         // Use blocked code

         if ( UPPER ) {

            // Compute inverse of upper triangular matrix

            for (J = 1; NB < 0 ? J >= N : J <= N; J += NB) { // 20
               JB = min( NB, N-J+1 );

               // Compute rows 1:j-1 of current block column

               dtrmm('Left', 'Upper', 'No transpose', DIAG, J-1, JB, ONE, A, LDA, A( 1, J ), LDA );
               dtrsm('Right', 'Upper', 'No transpose', DIAG, J-1, JB, -ONE, A( J, J ), LDA, A( 1, J ), LDA );

               // Compute inverse of current diagonal block

               dtrti2('Upper', DIAG, JB, A( J, J ), LDA, INFO );
            } // 20
         } else {

            // Compute inverse of lower triangular matrix

            NN = ( ( N-1 ) / NB )*NB + 1;
            for (J = NN; -NB < 0 ? J >= 1 : J <= 1; J += -NB) { // 30
               JB = min( NB, N-J+1 );
               if ( J+JB <= N ) {

                  // Compute rows j+jb:n of current block column

                  dtrmm('Left', 'Lower', 'No transpose', DIAG, N-J-JB+1, JB, ONE, A( J+JB, J+JB ), LDA, A( J+JB, J ), LDA );
                  dtrsm('Right', 'Lower', 'No transpose', DIAG, N-J-JB+1, JB, -ONE, A( J, J ), LDA, A( J+JB, J ), LDA );
               }

               // Compute inverse of current diagonal block

               dtrti2('Lower', DIAG, JB, A( J, J ), LDA, INFO );
            } // 30
         }
      }

      return;
      }
