      SUBROUTINE SLAPTM( N, NRHS, ALPHA, D, E, X, LDX, BETA, B, LDB )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDB, LDX, N, NRHS;
      REAL               ALPHA, BETA
      // ..
      // .. Array Arguments ..
      REAL               B( LDB, * ), D( * ), E( * ), X( LDX, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      // ..
      // .. Executable Statements ..

      IF( N.EQ.0 ) RETURN

      // Multiply B by BETA if BETA.NE.1.

      if ( BETA.EQ.ZERO ) {
         DO 20 J = 1, NRHS
            DO 10 I = 1, N
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      } else if ( BETA.EQ.-ONE ) {
         DO 40 J = 1, NRHS
            DO 30 I = 1, N
               B( I, J ) = -B( I, J )
   30       CONTINUE
   40    CONTINUE
      }

      if ( ALPHA.EQ.ONE ) {

         // Compute B := B + A*X

         DO 60 J = 1, NRHS
            if ( N.EQ.1 ) {
               B( 1, J ) = B( 1, J ) + D( 1 )*X( 1, J )
            } else {
               B( 1, J ) = B( 1, J ) + D( 1 )*X( 1, J ) + E( 1 )*X( 2, J )                B( N, J ) = B( N, J ) + E( N-1 )*X( N-1, J ) + D( N )*X( N, J )
               DO 50 I = 2, N - 1
                  B( I, J ) = B( I, J ) + E( I-1 )*X( I-1, J ) + D( I )*X( I, J ) + E( I )*X( I+1, J )
   50          CONTINUE
            }
   60    CONTINUE
      } else if ( ALPHA.EQ.-ONE ) {

         // Compute B := B - A*X

         DO 80 J = 1, NRHS
            if ( N.EQ.1 ) {
               B( 1, J ) = B( 1, J ) - D( 1 )*X( 1, J )
            } else {
               B( 1, J ) = B( 1, J ) - D( 1 )*X( 1, J ) - E( 1 )*X( 2, J )                B( N, J ) = B( N, J ) - E( N-1 )*X( N-1, J ) - D( N )*X( N, J )
               DO 70 I = 2, N - 1
                  B( I, J ) = B( I, J ) - E( I-1 )*X( I-1, J ) - D( I )*X( I, J ) - E( I )*X( I+1, J )
   70          CONTINUE
            }
   80    CONTINUE
      }
      RETURN

      // End of SLAPTM

      }
