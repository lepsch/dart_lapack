      SUBROUTINE CTRT03( UPLO, TRANS, DIAG, N, NRHS, A, LDA, SCALE, CNORM, TSCAL, X, LDX, B, LDB, WORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, TRANS, UPLO;
      int                LDA, LDB, LDX, N, NRHS;
      REAL               RESID, SCALE, TSCAL
      // ..
      // .. Array Arguments ..
      REAL               CNORM( * )
      COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * ), X( LDX, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                IX, J;
      REAL               EPS, ERR, SMLNUM, TNORM, XNORM, XSCAL
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ICAMAX;
      REAL               SLAMCH
      // EXTERNAL LSAME, ICAMAX, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CCOPY, CSSCAL, CTRMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CMPLX, MAX, REAL
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0

      IF( N.LE.0 .OR. NRHS.LE.0 ) THEN
         RESID = ZERO
         RETURN
      END IF
      EPS = SLAMCH( 'Epsilon' )
      SMLNUM = SLAMCH( 'Safe minimum' )

      // Compute the norm of the triangular matrix A using the column
      // norms already computed by CLATRS.

      TNORM = ZERO
      IF( LSAME( DIAG, 'N' ) ) THEN
         DO 10 J = 1, N
            TNORM = MAX( TNORM, TSCAL*ABS( A( J, J ) )+CNORM( J ) )
   10    CONTINUE
      ELSE
         DO 20 J = 1, N
            TNORM = MAX( TNORM, TSCAL+CNORM( J ) )
   20    CONTINUE
      END IF

      // Compute the maximum over the number of right hand sides of
         // norm(op(A)*x - s*b) / ( norm(op(A)) * norm(x) * EPS ).

      RESID = ZERO
      DO 30 J = 1, NRHS
         CALL CCOPY( N, X( 1, J ), 1, WORK, 1 )
         IX = ICAMAX( N, WORK, 1 )
         XNORM = MAX( ONE, ABS( X( IX, J ) ) )
         XSCAL = ( ONE / XNORM ) / REAL( N )
         CALL CSSCAL( N, XSCAL, WORK, 1 )
         CALL CTRMV( UPLO, TRANS, DIAG, N, A, LDA, WORK, 1 )
         CALL CAXPY( N, CMPLX( -SCALE*XSCAL ), B( 1, J ), 1, WORK, 1 )
         IX = ICAMAX( N, WORK, 1 )
         ERR = TSCAL*ABS( WORK( IX ) )
         IX = ICAMAX( N, X( 1, J ), 1 )
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
   30 CONTINUE

      RETURN

      // End of CTRT03

      }
