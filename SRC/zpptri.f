      SUBROUTINE ZPPTRI( UPLO, N, AP, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         AP( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE;
      PARAMETER          ( ONE = 1.0D+0 )
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                J, JC, JJ, JJN;
      double             AJJ;
      // ..
      // .. External Functions ..
      bool               LSAME;
      COMPLEX*16         ZDOTC
      // EXTERNAL LSAME, ZDOTC
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZDSCAL, ZHPR, ZTPMV, ZTPTRI
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZPPTRI', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // Invert the triangular Cholesky factor U or L.

      CALL ZTPTRI( UPLO, 'Non-unit', N, AP, INFO )
      IF( INFO.GT.0 ) RETURN
      IF( UPPER ) THEN

         // Compute the product inv(U) * inv(U)**H.

         JJ = 0
         DO 10 J = 1, N
            JC = JJ + 1
            JJ = JJ + J
            IF( J.GT.1 ) CALL ZHPR( 'Upper', J-1, ONE, AP( JC ), 1, AP )
            AJJ = DBLE( AP( JJ ) )
            CALL ZDSCAL( J, AJJ, AP( JC ), 1 )
   10    CONTINUE

      ELSE

         // Compute the product inv(L)**H * inv(L).

         JJ = 1
         DO 20 J = 1, N
            JJN = JJ + N - J + 1
            AP( JJ ) = DBLE( ZDOTC( N-J+1, AP( JJ ), 1, AP( JJ ), 1 ) )
            IF( J.LT.N ) CALL ZTPMV( 'Lower', 'Conjugate transpose', 'Non-unit', N-J, AP( JJN ), AP( JJ+1 ), 1 )
            JJ = JJN
   20    CONTINUE
      END IF

      RETURN

      // End of ZPPTRI

      END
