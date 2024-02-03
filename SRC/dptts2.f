      SUBROUTINE DPTTS2( N, NRHS, D, E, B, LDB )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      double             B( LDB, * ), D( * ), E( * );
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      int                I, J
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSCAL
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.LE.1 ) THEN
         IF( N.EQ.1 ) CALL DSCAL( NRHS, 1.D0 / D( 1 ), B, LDB )
         RETURN
      END IF
*
*     Solve A * X = B using the factorization A = L*D*L**T,
*     overwriting each right hand side vector with its solution.
*
      DO 30 J = 1, NRHS
*
*           Solve L * x = b.
*
         DO 10 I = 2, N
            B( I, J ) = B( I, J ) - B( I-1, J )*E( I-1 )
   10    CONTINUE
*
*           Solve D * L**T * x = b.
*
         B( N, J ) = B( N, J ) / D( N )
         DO 20 I = N - 1, 1, -1
            B( I, J ) = B( I, J ) / D( I ) - B( I+1, J )*E( I )
   20    CONTINUE
   30 CONTINUE
*
      RETURN
*
*     End of DPTTS2
*
      END
