      SUBROUTINE ZGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, KL, KU, LDAB, M, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX*16         AB( LDAB, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         ONE, ZERO
      const              ONE = ( 1.0D+0, 0.0D+0 ), ZERO = ( 0.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, J, JP, JU, KM, KV;
      // ..
      // .. External Functions ..
      int                IZAMAX;
      // EXTERNAL IZAMAX
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZGERU, ZSCAL, ZSWAP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // KV is the number of superdiagonals in the factor U, allowing for
      // fill-in.

      KV = KU + KL

      // Test the input parameters.

      INFO = 0
      if ( M.LT.0 ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( KL.LT.0 ) {
         INFO = -3
      } else if ( KU.LT.0 ) {
         INFO = -4
      } else if ( LDAB.LT.KL+KV+1 ) {
         INFO = -6
      }
      if ( INFO.NE.0 ) {
         xerbla('ZGBTF2', -INFO );
         RETURN
      }

      // Quick return if possible

      if (M.EQ.0 .OR. N.EQ.0) RETURN;

      // Gaussian elimination with partial pivoting

      // Set fill-in elements in columns KU+2 to KV to zero.

      DO 20 J = KU + 2, MIN( KV, N )
         for (I = KV - J + 2; I <= KL; I++) { // 10
            AB( I, J ) = ZERO
         } // 10
      } // 20

      // JU is the index of the last column affected by the current stage
      // of the factorization.

      JU = 1

      DO 40 J = 1, MIN( M, N )

         // Set fill-in elements in column J+KV to zero.

         if ( J+KV.LE.N ) {
            for (I = 1; I <= KL; I++) { // 30
               AB( I, J+KV ) = ZERO
            } // 30
         }

         // Find pivot and test for singularity. KM is the number of
         // subdiagonal elements in the current column.

         KM = MIN( KL, M-J )
         JP = IZAMAX( KM+1, AB( KV+1, J ), 1 )
         IPIV( J ) = JP + J - 1
         if ( AB( KV+JP, J ).NE.ZERO ) {
            JU = MAX( JU, MIN( J+KU+JP-1, N ) )

            // Apply interchange to columns J to JU.

            if (JP.NE.1) CALL ZSWAP( JU-J+1, AB( KV+JP, J ), LDAB-1, AB( KV+1, J ), LDAB-1 );
            if ( KM.GT.0 ) {

               // Compute multipliers.

               zscal(KM, ONE / AB( KV+1, J ), AB( KV+2, J ), 1 );

               // Update trailing submatrix within the band.

               if (JU.GT.J) CALL ZGERU( KM, JU-J, -ONE, AB( KV+2, J ), 1, AB( KV, J+1 ), LDAB-1, AB( KV+1, J+1 ), LDAB-1 );
            }
         } else {

            // If pivot is zero, set INFO to the index of the pivot
            // unless a zero pivot has already been found.

            if (INFO.EQ.0) INFO = J;
         }
      } // 40
      RETURN

      // End of ZGBTF2

      }
