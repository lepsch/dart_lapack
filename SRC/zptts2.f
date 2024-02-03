      SUBROUTINE ZPTTS2( IUPLO, N, NRHS, D, E, B, LDB )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                IUPLO, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      double             D( * );
      COMPLEX*16         B( LDB, * ), E( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      int                I, J
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZDSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DCONJG
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.LE.1 ) THEN
         IF( N.EQ.1 ) CALL ZDSCAL( NRHS, 1.D0 / D( 1 ), B, LDB )
         RETURN
      END IF
*
      IF( IUPLO.EQ.1 ) THEN
*
*        Solve A * X = B using the factorization A = U**H *D*U,
*        overwriting each right hand side vector with its solution.
*
         IF( NRHS.LE.2 ) THEN
            J = 1
   10       CONTINUE
*
*           Solve U**H * x = b.
*
            DO 20 I = 2, N
               B( I, J ) = B( I, J ) - B( I-1, J )*DCONJG( E( I-1 ) )
   20       CONTINUE
*
*           Solve D * U * x = b.
*
            DO 30 I = 1, N
               B( I, J ) = B( I, J ) / D( I )
   30       CONTINUE
            DO 40 I = N - 1, 1, -1
               B( I, J ) = B( I, J ) - B( I+1, J )*E( I )
   40       CONTINUE
            IF( J.LT.NRHS ) THEN
               J = J + 1
               GO TO 10
            END IF
         ELSE
            DO 70 J = 1, NRHS
*
*              Solve U**H * x = b.
*
               DO 50 I = 2, N
                  B( I, J ) = B( I, J ) - B( I-1, J )*DCONJG( E( I-1 ) )
   50          CONTINUE
*
*              Solve D * U * x = b.
*
               B( N, J ) = B( N, J ) / D( N )
               DO 60 I = N - 1, 1, -1
                  B( I, J ) = B( I, J ) / D( I ) - B( I+1, J )*E( I )
   60          CONTINUE
   70       CONTINUE
         END IF
      ELSE
*
*        Solve A * X = B using the factorization A = L*D*L**H,
*        overwriting each right hand side vector with its solution.
*
         IF( NRHS.LE.2 ) THEN
            J = 1
   80       CONTINUE
*
*           Solve L * x = b.
*
            DO 90 I = 2, N
               B( I, J ) = B( I, J ) - B( I-1, J )*E( I-1 )
   90       CONTINUE
*
*           Solve D * L**H * x = b.
*
            DO 100 I = 1, N
               B( I, J ) = B( I, J ) / D( I )
  100       CONTINUE
            DO 110 I = N - 1, 1, -1
               B( I, J ) = B( I, J ) - B( I+1, J )*DCONJG( E( I ) )
  110       CONTINUE
            IF( J.LT.NRHS ) THEN
               J = J + 1
               GO TO 80
            END IF
         ELSE
            DO 140 J = 1, NRHS
*
*              Solve L * x = b.
*
               DO 120 I = 2, N
                  B( I, J ) = B( I, J ) - B( I-1, J )*E( I-1 )
  120          CONTINUE
*
*              Solve D * L**H * x = b.
*
               B( N, J ) = B( N, J ) / D( N )
               DO 130 I = N - 1, 1, -1
                  B( I, J ) = B( I, J ) / D( I ) - B( I+1, J )*DCONJG( E( I ) )
  130          CONTINUE
  140       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of ZPTTS2
*
      END
