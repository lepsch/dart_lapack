      SUBROUTINE SBDT04( UPLO, N, D, E, S, NS, U, LDU, VT, LDVT, WORK, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String             UPLO;
      int                LDU, LDVT, N, NS;
      REAL               RESID
      // ..
      // .. Array Arguments ..
      REAL               D( * ), E( * ), S( * ), U( LDU, * ), VT( LDVT, * ), WORK( * )
      // ..
*
* ======================================================================
*
      // .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      // ..
      // .. Local Scalars ..
      int                I, J, K;
      REAL               BNORM, EPS
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ISAMAX;
      REAL               SASUM, SLAMCH
      // EXTERNAL LSAME, ISAMAX, SASUM, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, REAL, MAX, MIN
      // ..
      // .. Executable Statements ..
*
      // Quick return if possible.
*
      RESID = ZERO
      IF( N.LE.0 .OR. NS.LE.0 ) RETURN
*
      EPS = SLAMCH( 'Precision' )
*
      // Compute S - U' * B * V.
*
      BNORM = ZERO
*
      IF( LSAME( UPLO, 'U' ) ) THEN
*
         // B is upper bidiagonal.
*
         K = 0
         DO 20 I = 1, NS
            DO 10 J = 1, N-1
               K = K + 1
               WORK( K ) = D( J )*VT( I, J ) + E( J )*VT( I, J+1 )
   10       CONTINUE
            K = K + 1
            WORK( K ) = D( N )*VT( I, N )
   20    CONTINUE
         BNORM = ABS( D( 1 ) )
         DO 30 I = 2, N
            BNORM = MAX( BNORM, ABS( D( I ) )+ABS( E( I-1 ) ) )
   30    CONTINUE
      ELSE
*
         // B is lower bidiagonal.
*
         K = 0
         DO 50 I = 1, NS
            K = K + 1
            WORK( K ) = D( 1 )*VT( I, 1 )
            DO 40 J = 1, N-1
               K = K + 1
               WORK( K ) = E( J )*VT( I, J ) + D( J+1 )*VT( I, J+1 )
   40       CONTINUE
   50    CONTINUE
         BNORM = ABS( D( N ) )
         DO 60 I = 1, N-1
            BNORM = MAX( BNORM, ABS( D( I ) )+ABS( E( I ) ) )
   60    CONTINUE
      END IF
*
      CALL SGEMM( 'T', 'N', NS, NS, N, -ONE, U, LDU, WORK( 1 ), N, ZERO, WORK( 1+N*NS ), NS )
*
      // norm(S - U' * B * V)
*
      K = N*NS
      DO 70 I = 1, NS
         WORK( K+I ) =  WORK( K+I ) + S( I )
         RESID = MAX( RESID, SASUM( NS, WORK( K+1 ), 1 ) )
         K = K + NS
   70 CONTINUE
*
      IF( BNORM.LE.ZERO ) THEN
         IF( RESID.NE.ZERO ) RESID = ONE / EPS
      ELSE
         IF( BNORM.GE.RESID ) THEN
            RESID = ( RESID / BNORM ) / ( REAL( N )*EPS )
         ELSE
            IF( BNORM.LT.ONE ) THEN
               RESID = ( MIN( RESID, REAL( N )*BNORM ) / BNORM ) / ( REAL( N )*EPS )
            ELSE
               RESID = MIN( RESID / BNORM, REAL( N ) ) / ( REAL( N )*EPS )
            END IF
         END IF
      END IF
*
      RETURN
*
      // End of SBDT04
*
      END
