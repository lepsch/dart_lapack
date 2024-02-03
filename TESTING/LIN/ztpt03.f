      SUBROUTINE ZTPT03( UPLO, TRANS, DIAG, N, NRHS, AP, SCALE, CNORM, TSCAL, X, LDX, B, LDB, WORK, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                LDB, LDX, N, NRHS
      double             RESID, SCALE, TSCAL;
*     ..
*     .. Array Arguments ..
      double             CNORM( * );
      COMPLEX*16         AP( * ), B( LDB, * ), WORK( * ), X( LDX, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ONE, ZERO;
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      int                IX, J, JJ
      double             EPS, ERR, SMLNUM, TNORM, XNORM, XSCAL;
*     ..
*     .. External Functions ..
      bool               LSAME;
      int                IZAMAX
      double             DLAMCH;
      EXTERNAL           LSAME, IZAMAX, DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZAXPY, ZCOPY, ZDSCAL, ZTPMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, MAX
*     ..
*     .. Executable Statements ..
*
*     Quick exit if N = 0.
*
      IF( N.LE.0 .OR. NRHS.LE.0 ) THEN
         RESID = ZERO
         RETURN
      END IF
      EPS = DLAMCH( 'Epsilon' )
      SMLNUM = DLAMCH( 'Safe minimum' )
*
*     Compute the norm of the triangular matrix A using the column
*     norms already computed by ZLATPS.
*
      TNORM = 0.D0
      IF( LSAME( DIAG, 'N' ) ) THEN
         IF( LSAME( UPLO, 'U' ) ) THEN
            JJ = 1
            DO 10 J = 1, N
               TNORM = MAX( TNORM, TSCAL*ABS( AP( JJ ) )+CNORM( J ) )
               JJ = JJ + J
   10       CONTINUE
         ELSE
            JJ = 1
            DO 20 J = 1, N
               TNORM = MAX( TNORM, TSCAL*ABS( AP( JJ ) )+CNORM( J ) )
               JJ = JJ + N - J + 1
   20       CONTINUE
         END IF
      ELSE
         DO 30 J = 1, N
            TNORM = MAX( TNORM, TSCAL+CNORM( J ) )
   30    CONTINUE
      END IF
*
*     Compute the maximum over the number of right hand sides of
*        norm(op(A)*x - s*b) / ( norm(A) * norm(x) * EPS ).
*
      RESID = ZERO
      DO 40 J = 1, NRHS
         CALL ZCOPY( N, X( 1, J ), 1, WORK, 1 )
         IX = IZAMAX( N, WORK, 1 )
         XNORM = MAX( ONE, ABS( X( IX, J ) ) )
         XSCAL = ( ONE / XNORM ) / DBLE( N )
         CALL ZDSCAL( N, XSCAL, WORK, 1 )
         CALL ZTPMV( UPLO, TRANS, DIAG, N, AP, WORK, 1 )
         CALL ZAXPY( N, DCMPLX( -SCALE*XSCAL ), B( 1, J ), 1, WORK, 1 )
         IX = IZAMAX( N, WORK, 1 )
         ERR = TSCAL*ABS( WORK( IX ) )
         IX = IZAMAX( N, X( 1, J ), 1 )
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
*     End of ZTPT03
*
      END
