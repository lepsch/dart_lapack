      SUBROUTINE STPLQT( M, N, L, MB, A, LDA, B, LDB, T, LDT, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     INFO, LDA, LDB, LDT, N, M, L, MB;
      // ..
      // .. Array Arguments ..
      REAL    A( LDA, * ), B( LDB, * ), T( LDT, * ), WORK( * )
      // ..

* =====================================================================

      // ..
      // .. Local Scalars ..
      int        I, IB, LB, NB, IINFO;
      // ..
      // .. External Subroutines ..
      // EXTERNAL STPLQT2, STPRFB, XERBLA
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      if ( M.LT.0 ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( L.LT.0 .OR. (L.GT.MIN(M,N) .AND. MIN(M,N).GE.0)) {
         INFO = -3
      } else if ( MB.LT.1 .OR. (MB.GT.M .AND. M.GT.0)) {
         INFO = -4
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -6
      } else if ( LDB.LT.MAX( 1, M ) ) {
         INFO = -8
      } else if ( LDT.LT.MB ) {
         INFO = -10
      }
      if ( INFO.NE.0 ) {
         xerbla('STPLQT', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( M.EQ.0 .OR. N.EQ.0 ) RETURN

      DO I = 1, M, MB

      // Compute the QR factorization of the current block

         IB = MIN( M-I+1, MB )
         NB = MIN( N-L+I+IB-1, N )
         if ( I.GE.L ) {
            LB = 0
         } else {
            LB = NB-N+L-I+1
         }

         stplqt2(IB, NB, LB, A(I,I), LDA, B( I, 1 ), LDB, T(1, I ), LDT, IINFO );

      // Update by applying H**T to B(I+IB:M,:) from the right

         if ( I+IB.LE.M ) {
            stprfb('R', 'N', 'F', 'R', M-I-IB+1, NB, IB, LB, B( I, 1 ), LDB, T( 1, I ), LDT, A( I+IB, I ), LDA, B( I+IB, 1 ), LDB, WORK, M-I-IB+1);
         }
      END DO
      RETURN

      // End of STPLQT

      }
