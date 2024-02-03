      SUBROUTINE DGBT02( TRANS, M, N, KL, KU, NRHS, A, LDA, X, LDX, B, LDB, RWORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                KL, KU, LDA, LDB, LDX, M, N, NRHS;
      double             RESID;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), B( LDB, * ), X( LDX, * ), RWORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                I1, I2, J, KD, N1;
      double             ANORM, BNORM, EPS, TEMP, XNORM;
      // ..
      // .. External Functions ..
      bool               DISNAN, LSAME;
      double             DASUM, DLAMCH;
      // EXTERNAL DASUM, DISNAN, DLAMCH, LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGBMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Quick return if N = 0 pr NRHS = 0

      IF( M.LE.0 .OR. N.LE.0 .OR. NRHS.LE.0 ) THEN
         RESID = ZERO
         RETURN
      END IF

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = DLAMCH( 'Epsilon' )
      ANORM = ZERO
      IF( LSAME( TRANS, 'N' ) ) THEN

         // Find norm1(A).

         KD = KU + 1
         DO 10 J = 1, N
            I1 = MAX( KD+1-J, 1 )
            I2 = MIN( KD+M-J, KL+KD )
            IF( I2.GE.I1 ) THEN
               TEMP = DASUM( I2-I1+1, A( I1, J ), 1 )
               IF( ANORM.LT.TEMP .OR. DISNAN( TEMP ) ) ANORM = TEMP
            END IF
   10    CONTINUE
      ELSE

         // Find normI(A).

         DO 12 I1 = 1, M
            RWORK( I1 ) = ZERO
   12    CONTINUE
         DO 16 J = 1, N
            KD = KU + 1 - J
            DO 14 I1 = MAX( 1, J-KU ), MIN( M, J+KL )
               RWORK( I1 ) = RWORK( I1 ) + ABS( A( KD+I1, J ) )
   14       CONTINUE
   16    CONTINUE
         DO 18 I1 = 1, M
            TEMP = RWORK( I1 )
            IF( ANORM.LT.TEMP .OR. DISNAN( TEMP ) ) ANORM = TEMP
   18    CONTINUE
      END IF
      IF( ANORM.LE.ZERO ) THEN
         RESID = ONE / EPS
         RETURN
      END IF

      IF( LSAME( TRANS, 'T' ) .OR. LSAME( TRANS, 'C' ) ) THEN
         N1 = N
      ELSE
         N1 = M
      END IF

      // Compute B - op(A)*X

      DO 20 J = 1, NRHS
         CALL DGBMV( TRANS, M, N, KL, KU, -ONE, A, LDA, X( 1, J ), 1, ONE, B( 1, J ), 1 )
   20 CONTINUE

      // Compute the maximum over the number of right hand sides of
         // norm(B - op(A)*X) / ( norm(op(A)) * norm(X) * EPS ).

      RESID = ZERO
      DO 30 J = 1, NRHS
         BNORM = DASUM( N1, B( 1, J ), 1 )
         XNORM = DASUM( N1, X( 1, J ), 1 )
         IF( XNORM.LE.ZERO ) THEN
            RESID = ONE / EPS
         ELSE
            RESID = MAX( RESID, ( ( BNORM / ANORM ) / XNORM ) / EPS )
         END IF
   30 CONTINUE

      RETURN

      // End of DGBT02

      }
