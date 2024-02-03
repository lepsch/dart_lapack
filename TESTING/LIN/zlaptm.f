      SUBROUTINE ZLAPTM( UPLO, N, NRHS, ALPHA, D, E, X, LDX, BETA, B, LDB )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String             UPLO;
      int                LDB, LDX, N, NRHS;
      double             ALPHA, BETA;
      // ..
      // .. Array Arguments ..
      double             D( * );
      COMPLEX*16         B( LDB, * ), E( * ), X( LDX, * )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      double             ONE, ZERO;
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      // ..
      // .. Local Scalars ..
      int                I, J;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCONJG
      // ..
      // .. Executable Statements ..
*
      IF( N.EQ.0 ) RETURN
*
      IF( BETA.EQ.ZERO ) THEN
         DO 20 J = 1, NRHS
            DO 10 I = 1, N
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE IF( BETA.EQ.-ONE ) THEN
         DO 40 J = 1, NRHS
            DO 30 I = 1, N
               B( I, J ) = -B( I, J )
   30       CONTINUE
   40    CONTINUE
      END IF
*
      IF( ALPHA.EQ.ONE ) THEN
         IF( LSAME( UPLO, 'U' ) ) THEN
*
            // Compute B := B + A*X, where E is the superdiagonal of A.
*
            DO 60 J = 1, NRHS
               IF( N.EQ.1 ) THEN
                  B( 1, J ) = B( 1, J ) + D( 1 )*X( 1, J )
               ELSE
                  B( 1, J ) = B( 1, J ) + D( 1 )*X( 1, J ) + E( 1 )*X( 2, J )                   B( N, J ) = B( N, J ) + DCONJG( E( N-1 ) )* X( N-1, J ) + D( N )*X( N, J )
                  DO 50 I = 2, N - 1
                     B( I, J ) = B( I, J ) + DCONJG( E( I-1 ) )* X( I-1, J ) + D( I )*X( I, J ) + E( I )*X( I+1, J )
   50             CONTINUE
               END IF
   60       CONTINUE
         ELSE
*
            // Compute B := B + A*X, where E is the subdiagonal of A.
*
            DO 80 J = 1, NRHS
               IF( N.EQ.1 ) THEN
                  B( 1, J ) = B( 1, J ) + D( 1 )*X( 1, J )
               ELSE
                  B( 1, J ) = B( 1, J ) + D( 1 )*X( 1, J ) + DCONJG( E( 1 ) )*X( 2, J )                   B( N, J ) = B( N, J ) + E( N-1 )*X( N-1, J ) + D( N )*X( N, J )
                  DO 70 I = 2, N - 1
                     B( I, J ) = B( I, J ) + E( I-1 )*X( I-1, J ) + D( I )*X( I, J ) + DCONJG( E( I ) )*X( I+1, J )
   70             CONTINUE
               END IF
   80       CONTINUE
         END IF
      ELSE IF( ALPHA.EQ.-ONE ) THEN
         IF( LSAME( UPLO, 'U' ) ) THEN
*
            // Compute B := B - A*X, where E is the superdiagonal of A.
*
            DO 100 J = 1, NRHS
               IF( N.EQ.1 ) THEN
                  B( 1, J ) = B( 1, J ) - D( 1 )*X( 1, J )
               ELSE
                  B( 1, J ) = B( 1, J ) - D( 1 )*X( 1, J ) - E( 1 )*X( 2, J )                   B( N, J ) = B( N, J ) - DCONJG( E( N-1 ) )* X( N-1, J ) - D( N )*X( N, J )
                  DO 90 I = 2, N - 1
                     B( I, J ) = B( I, J ) - DCONJG( E( I-1 ) )* X( I-1, J ) - D( I )*X( I, J ) - E( I )*X( I+1, J )
   90             CONTINUE
               END IF
  100       CONTINUE
         ELSE
*
            // Compute B := B - A*X, where E is the subdiagonal of A.
*
            DO 120 J = 1, NRHS
               IF( N.EQ.1 ) THEN
                  B( 1, J ) = B( 1, J ) - D( 1 )*X( 1, J )
               ELSE
                  B( 1, J ) = B( 1, J ) - D( 1 )*X( 1, J ) - DCONJG( E( 1 ) )*X( 2, J )                   B( N, J ) = B( N, J ) - E( N-1 )*X( N-1, J ) - D( N )*X( N, J )
                  DO 110 I = 2, N - 1
                     B( I, J ) = B( I, J ) - E( I-1 )*X( I-1, J ) - D( I )*X( I, J ) - DCONJG( E( I ) )*X( I+1, J )
  110             CONTINUE
               END IF
  120       CONTINUE
         END IF
      END IF
      RETURN
*
      // End of ZLAPTM
*
      END
