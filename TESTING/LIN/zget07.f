      SUBROUTINE ZGET07( TRANS, N, NRHS, A, LDA, B, LDB, X, LDX, XACT, LDXACT, FERR, CHKFERR, BERR, RESLTS )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             TRANS;
      bool               CHKFERR;
      int                LDA, LDB, LDX, LDXACT, N, NRHS
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   BERR( * ), FERR( * ), RESLTS( * )
      COMPLEX*16         A( LDA, * ), B( LDB, * ), X( LDX, * ), XACT( LDXACT, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      bool               NOTRAN;
      int                I, IMAX, J, K
      DOUBLE PRECISION   AXBI, DIFF, EPS, ERRBND, OVFL, TMP, UNFL, XNORM
      COMPLEX*16         ZDUM
*     ..
*     .. External Functions ..
      bool               LSAME;
      int                IZAMAX
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, IZAMAX, DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DIMAG, MAX, MIN
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )
*     ..
*     .. Executable Statements ..
*
*     Quick exit if N = 0 or NRHS = 0.
*
      IF( N.LE.0 .OR. NRHS.LE.0 ) THEN
         RESLTS( 1 ) = ZERO
         RESLTS( 2 ) = ZERO
         RETURN
      END IF
*
      EPS = DLAMCH( 'Epsilon' )
      UNFL = DLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      NOTRAN = LSAME( TRANS, 'N' )
*
*     Test 1:  Compute the maximum of
*        norm(X - XACT) / ( norm(X) * FERR )
*     over all the vectors X and XACT using the infinity-norm.
*
      ERRBND = ZERO
      IF( CHKFERR ) THEN
         DO 30 J = 1, NRHS
            IMAX = IZAMAX( N, X( 1, J ), 1 )
            XNORM = MAX( CABS1( X( IMAX, J ) ), UNFL )
            DIFF = ZERO
            DO 10 I = 1, N
               DIFF = MAX( DIFF, CABS1( X( I, J )-XACT( I, J ) ) )
 10         CONTINUE
*
            IF( XNORM.GT.ONE ) THEN
               GO TO 20
            ELSE IF( DIFF.LE.OVFL*XNORM ) THEN
               GO TO 20
            ELSE
               ERRBND = ONE / EPS
               GO TO 30
            END IF
*
 20         CONTINUE
            IF( DIFF / XNORM.LE.FERR( J ) ) THEN
               ERRBND = MAX( ERRBND, ( DIFF / XNORM ) / FERR( J ) )
            ELSE
               ERRBND = ONE / EPS
            END IF
 30      CONTINUE
      END IF
      RESLTS( 1 ) = ERRBND
*
*     Test 2:  Compute the maximum of BERR / ( (n+1)*EPS + (*) ), where
*     (*) = (n+1)*UNFL / (min_i (abs(op(A))*abs(X) +abs(b))_i )
*
      DO 70 K = 1, NRHS
         DO 60 I = 1, N
            TMP = CABS1( B( I, K ) )
            IF( NOTRAN ) THEN
               DO 40 J = 1, N
                  TMP = TMP + CABS1( A( I, J ) )*CABS1( X( J, K ) )
   40          CONTINUE
            ELSE
               DO 50 J = 1, N
                  TMP = TMP + CABS1( A( J, I ) )*CABS1( X( J, K ) )
   50          CONTINUE
            END IF
            IF( I.EQ.1 ) THEN
               AXBI = TMP
            ELSE
               AXBI = MIN( AXBI, TMP )
            END IF
   60    CONTINUE
         TMP = BERR( K ) / ( ( N+1 )*EPS+( N+1 )*UNFL / MAX( AXBI, ( N+1 )*UNFL ) )
         IF( K.EQ.1 ) THEN
            RESLTS( 2 ) = TMP
         ELSE
            RESLTS( 2 ) = MAX( RESLTS( 2 ), TMP )
         END IF
   70 CONTINUE
*
      RETURN
*
*     End of ZGET07
*
      END
