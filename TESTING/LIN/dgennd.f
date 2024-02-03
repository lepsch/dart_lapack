      bool    FUNCTION DGENND (M, N, A, LDA);

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     M, N, LDA;
      // ..
      // .. Array Arguments ..
      double           A( LDA, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0D0 ;
      // ..
      // .. Local Scalars ..
      int     I, K;
      // ..
      // .. Intrinsics ..
      // INTRINSIC MIN
      // ..
      // .. Executable Statements ..
      K = MIN( M, N )
      DO I = 1, K
         if ( A( I, I ).LT.ZERO ) {
            DGENND = .FALSE.
            RETURN
         }
      END DO
      DGENND = .TRUE.
      RETURN
      }
