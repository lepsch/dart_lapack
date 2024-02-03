      SUBROUTINE DTBT03( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, SCALE, CNORM, TSCAL, X, LDX, B, LDB, WORK, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                KD, LDAB, LDB, LDX, N, NRHS
      DOUBLE PRECISION   RESID, SCALE, TSCAL
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   AB( LDAB, * ), B( LDB, * ), CNORM( * ), WORK( * ), X( LDX, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      int                IX, J
      DOUBLE PRECISION   BIGNUM, EPS, ERR, SMLNUM, TNORM, XNORM, XSCAL
*     ..
*     .. External Functions ..
      bool               LSAME;
      int                IDAMAX
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, IDAMAX, DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DAXPY, DCOPY, DSCAL, DTBMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX
*     ..
*     .. Executable Statements ..
*
*     Quick exit if N = 0
*
      IF( N.LE.0 .OR. NRHS.LE.0 ) THEN
         RESID = ZERO
         RETURN
      END IF
      EPS = DLAMCH( 'Epsilon' )
      SMLNUM = DLAMCH( 'Safe minimum' )
      BIGNUM = ONE / SMLNUM
*
*     Compute the norm of the triangular matrix A using the column
*     norms already computed by DLATBS.
*
      TNORM = ZERO
      IF( LSAME( DIAG, 'N' ) ) THEN
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 10 J = 1, N
               TNORM = MAX( TNORM, TSCAL*ABS( AB( KD+1, J ) )+ CNORM( J ) )
   10       CONTINUE
         ELSE
            DO 20 J = 1, N
               TNORM = MAX( TNORM, TSCAL*ABS( AB( 1, J ) )+CNORM( J ) )
   20       CONTINUE
         END IF
      ELSE
         DO 30 J = 1, N
            TNORM = MAX( TNORM, TSCAL+CNORM( J ) )
   30    CONTINUE
      END IF
*
*     Compute the maximum over the number of right hand sides of
*        norm(op(A)*x - s*b) / ( norm(op(A)) * norm(x) * EPS ).
*
      RESID = ZERO
      DO 40 J = 1, NRHS
         CALL DCOPY( N, X( 1, J ), 1, WORK, 1 )
         IX = IDAMAX( N, WORK, 1 )
         XNORM = MAX( ONE, ABS( X( IX, J ) ) )
         XSCAL = ( ONE / XNORM ) / DBLE( KD+1 )
         CALL DSCAL( N, XSCAL, WORK, 1 )
         CALL DTBMV( UPLO, TRANS, DIAG, N, KD, AB, LDAB, WORK, 1 )
         CALL DAXPY( N, -SCALE*XSCAL, B( 1, J ), 1, WORK, 1 )
         IX = IDAMAX( N, WORK, 1 )
         ERR = TSCAL*ABS( WORK( IX ) )
         IX = IDAMAX( N, X( 1, J ), 1 )
         XNORM = ABS( X( IX, J ) )
         IF( ERR*SMLNUM.LE.XNORM ) THEN
            IF( XNORM.GT.ZERO ) ERR = ERR / XNORM
         ELSE
            IF( ERR.GT.ZERO ) ERR = ONE / EPS
         END IF
         IF( ERR*SMLNUM.LE.TNORM ) THEN
            IF( TNORM.GT.ZERO ) ERR = ERR / TNORM
         ELSE
            IF( ERR.GT.ZERO ) ERR = ONE / EPS
         END IF
         RESID = MAX( RESID, ERR )
   40 CONTINUE
*
      RETURN
*
*     End of DTBT03
*
      END
