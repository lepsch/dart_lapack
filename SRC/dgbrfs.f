      SUBROUTINE DGBRFS( TRANS, N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB, IPIV, B, LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String             TRANS;
      int                INFO, KL, KU, LDAB, LDAFB, LDB, LDX, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * ), IWORK( * );
      double             AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ), BERR( * ), FERR( * ), WORK( * ), X( LDX, * );
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      int                ITMAX;
      PARAMETER          ( ITMAX = 5 )
      double             ZERO;
      PARAMETER          ( ZERO = 0.0D+0 )
      double             ONE;
      PARAMETER          ( ONE = 1.0D+0 )
      double             TWO;
      PARAMETER          ( TWO = 2.0D+0 )
      double             THREE;
      PARAMETER          ( THREE = 3.0D+0 )
      // ..
      // .. Local Scalars ..
      bool               NOTRAN;
      String             TRANST;
      int                COUNT, I, J, K, KASE, KK, NZ;
      double             EPS, LSTRES, S, SAFE1, SAFE2, SAFMIN, XK;
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DCOPY, DGBMV, DGBTRS, DLACN2, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH;
      // EXTERNAL LSAME, DLAMCH
      // ..
      // .. Executable Statements ..
*
      // Test the input parameters.
*
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDAB.LT.KL+KU+1 ) THEN
         INFO = -7
      ELSE IF( LDAFB.LT.2*KL+KU+1 ) THEN
         INFO = -9
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -12
      ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
         INFO = -14
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGBRFS', -INFO )
         RETURN
      END IF
*
      // Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 ) THEN
         DO 10 J = 1, NRHS
            FERR( J ) = ZERO
            BERR( J ) = ZERO
   10    CONTINUE
         RETURN
      END IF
*
      IF( NOTRAN ) THEN
         TRANST = 'T'
      ELSE
         TRANST = 'N'
      END IF
*
      // NZ = maximum number of nonzero elements in each row of A, plus 1
*
      NZ = MIN( KL+KU+2, N+1 )
      EPS = DLAMCH( 'Epsilon' )
      SAFMIN = DLAMCH( 'Safe minimum' )
      SAFE1 = NZ*SAFMIN
      SAFE2 = SAFE1 / EPS
*
      // Do for each right hand side
*
      DO 140 J = 1, NRHS
*
         COUNT = 1
         LSTRES = THREE
   20    CONTINUE
*
         // Loop until stopping criterion is satisfied.
*
         // Compute residual R = B - op(A) * X,
         // where op(A) = A, A**T, or A**H, depending on TRANS.
*
         CALL DCOPY( N, B( 1, J ), 1, WORK( N+1 ), 1 )
         CALL DGBMV( TRANS, N, N, KL, KU, -ONE, AB, LDAB, X( 1, J ), 1, ONE, WORK( N+1 ), 1 )
*
         // Compute componentwise relative backward error from formula
*
         // max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) )
*
         // where abs(Z) is the componentwise absolute value of the matrix
         // or vector Z.  If the i-th component of the denominator is less
        t // han SAFE2, then SAFE1 is added to the i-th components of the
         // numerator and denominator before dividing.
*
         DO 30 I = 1, N
            WORK( I ) = ABS( B( I, J ) )
   30    CONTINUE
*
         // Compute abs(op(A))*abs(X) + abs(B).
*
         IF( NOTRAN ) THEN
            DO 50 K = 1, N
               KK = KU + 1 - K
               XK = ABS( X( K, J ) )
               DO 40 I = MAX( 1, K-KU ), MIN( N, K+KL )
                  WORK( I ) = WORK( I ) + ABS( AB( KK+I, K ) )*XK
   40          CONTINUE
   50       CONTINUE
         ELSE
            DO 70 K = 1, N
               S = ZERO
               KK = KU + 1 - K
               DO 60 I = MAX( 1, K-KU ), MIN( N, K+KL )
                  S = S + ABS( AB( KK+I, K ) )*ABS( X( I, J ) )
   60          CONTINUE
               WORK( K ) = WORK( K ) + S
   70       CONTINUE
         END IF
         S = ZERO
         DO 80 I = 1, N
            IF( WORK( I ).GT.SAFE2 ) THEN
               S = MAX( S, ABS( WORK( N+I ) ) / WORK( I ) )
            ELSE
               S = MAX( S, ( ABS( WORK( N+I ) )+SAFE1 ) / ( WORK( I )+SAFE1 ) )
            END IF
   80    CONTINUE
         BERR( J ) = S
*
         // Test stopping criterion. Continue iterating if
            // 1) The residual BERR(J) is larger than machine epsilon, and
            // 2) BERR(J) decreased by at least a factor of 2 during the
               // last iteration, and
            // 3) At most ITMAX iterations tried.
*
         IF( BERR( J ).GT.EPS .AND. TWO*BERR( J ).LE.LSTRES .AND. COUNT.LE.ITMAX ) THEN
*
            // Update solution and try again.
*
            CALL DGBTRS( TRANS, N, KL, KU, 1, AFB, LDAFB, IPIV, WORK( N+1 ), N, INFO )
            CALL DAXPY( N, ONE, WORK( N+1 ), 1, X( 1, J ), 1 )
            LSTRES = BERR( J )
            COUNT = COUNT + 1
            GO TO 20
         END IF
*
         // Bound error from formula
*
         // norm(X - XTRUE) / norm(X) .le. FERR =
         // norm( abs(inv(op(A)))*
            // ( abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) ))) / norm(X)
*
         // where
           // norm(Z) is the magnitude of the largest component of Z
           // inv(op(A)) is the inverse of op(A)
           // abs(Z) is the componentwise absolute value of the matrix or
              // vector Z
           // NZ is the maximum number of nonzeros in any row of A, plus 1
           // EPS is machine epsilon
*
         // The i-th component of abs(R)+NZ*EPS*(abs(op(A))*abs(X)+abs(B))
         // is incremented by SAFE1 if the i-th component of
         // abs(op(A))*abs(X) + abs(B) is less than SAFE2.
*
         // Use DLACN2 to estimate the infinity-norm of the matrix
            // inv(op(A)) * diag(W),
         // where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) )))
*
         DO 90 I = 1, N
            IF( WORK( I ).GT.SAFE2 ) THEN
               WORK( I ) = ABS( WORK( N+I ) ) + NZ*EPS*WORK( I )
            ELSE
               WORK( I ) = ABS( WORK( N+I ) ) + NZ*EPS*WORK( I ) + SAFE1
            END IF
   90    CONTINUE
*
         KASE = 0
  100    CONTINUE
         CALL DLACN2( N, WORK( 2*N+1 ), WORK( N+1 ), IWORK, FERR( J ), KASE, ISAVE )
         IF( KASE.NE.0 ) THEN
            IF( KASE.EQ.1 ) THEN
*
               // Multiply by diag(W)*inv(op(A)**T).
*
               CALL DGBTRS( TRANST, N, KL, KU, 1, AFB, LDAFB, IPIV, WORK( N+1 ), N, INFO )
               DO 110 I = 1, N
                  WORK( N+I ) = WORK( N+I )*WORK( I )
  110          CONTINUE
            ELSE
*
               // Multiply by inv(op(A))*diag(W).
*
               DO 120 I = 1, N
                  WORK( N+I ) = WORK( N+I )*WORK( I )
  120          CONTINUE
               CALL DGBTRS( TRANS, N, KL, KU, 1, AFB, LDAFB, IPIV, WORK( N+1 ), N, INFO )
            END IF
            GO TO 100
         END IF
*
         // Normalize error.
*
         LSTRES = ZERO
         DO 130 I = 1, N
            LSTRES = MAX( LSTRES, ABS( X( I, J ) ) )
  130    CONTINUE
         IF( LSTRES.NE.ZERO ) FERR( J ) = FERR( J ) / LSTRES
*
  140 CONTINUE
*
      RETURN
*
      // End of DGBRFS
*
      END
