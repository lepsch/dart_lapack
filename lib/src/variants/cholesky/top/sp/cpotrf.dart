      import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';
import 'package:lapack/src/matrix.dart';

void cpotrf(final String UPLO, final int N, final Matrix<double> A_, final int LDA, final Box<int> INFO,) {
// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  final A = A_.dim(LDA);
      const              ONE = 1.0;
      bool               UPPER;
      int                J, JB, NB;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV;
      // EXTERNAL lsame, ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CHERK, CPOTRF2, CTRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

      // Test the input parameters.

      INFO.value = 0;
      UPPER = lsame( UPLO, 'U' );
      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO.value = -1;
      } else if ( N < 0 ) {
         INFO.value = -2;
      } else if ( LDA < max( 1, N ) ) {
         INFO.value = -4;
      }
      if ( INFO.value != 0 ) {
         xerbla('CPOTRF', -INFO.value );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Determine the block size for this environment.

      NB = ilaenv( 1, 'CPOTRF', UPLO, N, -1, -1, -1 );
      if ( NB <= 1 || NB >= N ) {

         // Use unblocked code.

         cpotrf2(UPLO, N, A, LDA, INFO.value );
      } else {

         // Use blocked code.

         if ( UPPER ) {

            // Compute the Cholesky factorization A = U'*U.

            for (J = 1; NB < 0 ? J >= N : J <= N; J += NB) { // 10

               JB = min( NB, N-J+1 );

               // Compute the current block.

               ctrsm('Left', 'Upper', 'Conjugate Transpose', 'Non-unit', J-1, JB, Complex.one, A( 1, 1 ), LDA, A( 1, J ), LDA );
                cherk('Upper', 'Conjugate Transpose', JB, J-1, -ONE, A( 1, J ), LDA, ONE, A( J, J ), LDA );

               // Update and factorize the current diagonal block and test
               // for non-positive-definiteness.

               cpotrf2('Upper', JB, A( J, J ), LDA, INFO.value );
               if (INFO.value != 0) GO TO 30;

            } // 10

         } else {

            // Compute the Cholesky factorization A = L*L'.

            for (J = 1; NB < 0 ? J >= N : J <= N; J += NB) { // 20

               JB = min( NB, N-J+1 );

               // Compute the current block.

               ctrsm('Right', 'Lower', 'Conjugate Transpose', 'Non-unit', JB, J-1, Complex.one, A( 1, 1 ), LDA, A( J, 1 ), LDA );
                cherk('Lower', 'No Transpose', JB, J-1, -ONE, A( J, 1 ), LDA, ONE, A( J, J ), LDA );

               // Update and factorize the current diagonal block and test
               // for non-positive-definiteness.

               cpotrf2('Lower', JB, A( J, J ), LDA, INFO.value );
               if (INFO.value != 0) GO TO 30;

            } // 20
         }
      }
      GO TO 40;

      } // 30
      INFO.value = INFO.value + J - 1;

      } // 40
      }
