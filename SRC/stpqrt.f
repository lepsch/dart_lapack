      SUBROUTINE STPQRT( M, N, L, NB, A, LDA, B, LDB, T, LDT, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INFO, LDA, LDB, LDT, N, M, L, NB;
      // ..
      // .. Array Arguments ..
      REAL A( LDA, * ), B( LDB, * ), T( LDT, * ), WORK( * )
      // ..

* =====================================================================

      // ..
      // .. Local Scalars ..
      int        I, IB, LB, MB, IINFO;
      // ..
      // .. External Subroutines ..
      // EXTERNAL STPQRT2, STPRFB, XERBLA
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      if ( M.LT.0 ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( L.LT.0 .OR. (L.GT.MIN(M,N) && MIN(M,N).GE.0)) {
         INFO = -3
      } else if ( NB.LT.1 .OR. (NB.GT.N && N.GT.0)) {
         INFO = -4
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -6
      } else if ( LDB.LT.MAX( 1, M ) ) {
         INFO = -8
      } else if ( LDT.LT.NB ) {
         INFO = -10
      }
      if ( INFO != 0 ) {
         xerbla('STPQRT', -INFO );
         RETURN
      }

      // Quick return if possible

      if (M == 0 .OR. N == 0) RETURN;

      DO I = 1, N, NB

      // Compute the QR factorization of the current block

         IB = MIN( N-I+1, NB )
         MB = MIN( M-L+I+IB-1, M )
         if ( I.GE.L ) {
            LB = 0
         } else {
            LB = MB-M+L-I+1
         }

         stpqrt2(MB, IB, LB, A(I,I), LDA, B( 1, I ), LDB, T(1, I ), LDT, IINFO );

      // Update by applying H^H to B(:,I+IB:N) from the left

         if ( I+IB.LE.N ) {
            stprfb('L', 'T', 'F', 'C', MB, N-I-IB+1, IB, LB, B( 1, I ), LDB, T( 1, I ), LDT, A( I, I+IB ), LDA, B( 1, I+IB ), LDB, WORK, IB );
         }
      }
      RETURN

      // End of STPQRT

      }
