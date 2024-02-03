      SUBROUTINE ZPST01( UPLO, N, A, LDA, AFAC, LDAFAC, PERM, LDPERM, PIV, RWORK, RESID, RANK )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double             RESID;
      int                LDA, LDAFAC, LDPERM, N, RANK;
      String             UPLO;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), AFAC( LDAFAC, * ), PERM( LDPERM, * )
      double             RWORK( * );
      int                PIV( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CZERO
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ) )
      // ..
      // .. Local Scalars ..
      COMPLEX*16         TC
      double             ANORM, EPS, TR;
      int                I, J, K;
      // ..
      // .. External Functions ..
      COMPLEX*16         ZDOTC
      double             DLAMCH, ZLANHE;
      bool               LSAME;
      // EXTERNAL ZDOTC, DLAMCH, ZLANHE, LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZHER, ZSCAL, ZTRMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCONJG, DIMAG
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0.

      IF( N.LE.0 ) THEN
         RESID = ZERO
         RETURN
      END IF

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = DLAMCH( 'Epsilon' )
      ANORM = ZLANHE( '1', UPLO, N, A, LDA, RWORK )
      IF( ANORM.LE.ZERO ) THEN
         RESID = ONE / EPS
         RETURN
      END IF

      // Check the imaginary parts of the diagonal elements and return with
      // an error code if any are nonzero.

      DO 100 J = 1, N
         IF( DIMAG( AFAC( J, J ) ).NE.ZERO ) THEN
            RESID = ONE / EPS
            RETURN
         END IF
  100 CONTINUE

      // Compute the product U'*U, overwriting U.

      IF( LSAME( UPLO, 'U' ) ) THEN

         IF( RANK.LT.N ) THEN
            DO 120 J = RANK + 1, N
               DO 110 I = RANK + 1, J
                  AFAC( I, J ) = CZERO
  110          CONTINUE
  120       CONTINUE
         END IF

         DO 130 K = N, 1, -1

            // Compute the (K,K) element of the result.

            TR = DBLE( ZDOTC( K, AFAC( 1, K ), 1, AFAC( 1, K ), 1 ) )
            AFAC( K, K ) = TR

            // Compute the rest of column K.

            CALL ZTRMV( 'Upper', 'Conjugate', 'Non-unit', K-1, AFAC, LDAFAC, AFAC( 1, K ), 1 )

  130    CONTINUE

      // Compute the product L*L', overwriting L.

      ELSE

         IF( RANK.LT.N ) THEN
            DO 150 J = RANK + 1, N
               DO 140 I = J, N
                  AFAC( I, J ) = CZERO
  140          CONTINUE
  150       CONTINUE
         END IF

         DO 160 K = N, 1, -1
            // Add a multiple of column K of the factor L to each of
            // columns K+1 through N.

            IF( K+1.LE.N ) CALL ZHER( 'Lower', N-K, ONE, AFAC( K+1, K ), 1, AFAC( K+1, K+1 ), LDAFAC )

            // Scale column K by the diagonal element.

            TC = AFAC( K, K )
            CALL ZSCAL( N-K+1, TC, AFAC( K, K ), 1 )
  160    CONTINUE

      END IF

         // Form P*L*L'*P' or P*U'*U*P'

      IF( LSAME( UPLO, 'U' ) ) THEN

         DO 180 J = 1, N
            DO 170 I = 1, N
               IF( PIV( I ).LE.PIV( J ) ) THEN
                  IF( I.LE.J ) THEN
                     PERM( PIV( I ), PIV( J ) ) = AFAC( I, J )
                  ELSE
                     PERM( PIV( I ), PIV( J ) ) = DCONJG( AFAC( J, I ) )
                  END IF
               END IF
  170       CONTINUE
  180    CONTINUE


      ELSE

         DO 200 J = 1, N
            DO 190 I = 1, N
               IF( PIV( I ).GE.PIV( J ) ) THEN
                  IF( I.GE.J ) THEN
                     PERM( PIV( I ), PIV( J ) ) = AFAC( I, J )
                  ELSE
                     PERM( PIV( I ), PIV( J ) ) = DCONJG( AFAC( J, I ) )
                  END IF
               END IF
  190       CONTINUE
  200    CONTINUE

      END IF

      // Compute the difference  P*L*L'*P' - A (or P*U'*U*P' - A).

      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 220 J = 1, N
            DO 210 I = 1, J - 1
               PERM( I, J ) = PERM( I, J ) - A( I, J )
  210       CONTINUE
            PERM( J, J ) = PERM( J, J ) - DBLE( A( J, J ) )
  220    CONTINUE
      ELSE
         DO 240 J = 1, N
            PERM( J, J ) = PERM( J, J ) - DBLE( A( J, J ) )
            DO 230 I = J + 1, N
               PERM( I, J ) = PERM( I, J ) - A( I, J )
  230       CONTINUE
  240    CONTINUE
      END IF

      // Compute norm( P*L*L'P - A ) / ( N * norm(A) * EPS ), or
      // ( P*U'*U*P' - A )/ ( N * norm(A) * EPS ).

      RESID = ZLANHE( '1', UPLO, N, PERM, LDAFAC, RWORK )

      RESID = ( ( RESID / DBLE( N ) ) / ANORM ) / EPS

      RETURN

      // End of ZPST01

      END
