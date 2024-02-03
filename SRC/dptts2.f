      SUBROUTINE DPTTS2( N, NRHS, D, E, B, LDB )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      double             B( LDB, * ), D( * ), E( * );
      // ..

*  =====================================================================

      // .. Local Scalars ..
      int                I, J;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSCAL
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if ( N.LE.1 ) {
         if (N == 1) CALL DSCAL( NRHS, 1.D0 / D( 1 ), B, LDB );
         RETURN
      }

      // Solve A * X = B using the factorization A = L*D*L**T,
      // overwriting each right hand side vector with its solution.

      for (J = 1; J <= NRHS; J++) { // 30

            // Solve L * x = b.

         for (I = 2; I <= N; I++) { // 10
            B( I, J ) = B( I, J ) - B( I-1, J )*E( I-1 )
         } // 10

            // Solve D * L**T * x = b.

         B( N, J ) = B( N, J ) / D( N )
         DO 20 I = N - 1, 1, -1
            B( I, J ) = B( I, J ) / D( I ) - B( I+1, J )*E( I )
         } // 20
      } // 30

      RETURN

      // End of DPTTS2

      }
