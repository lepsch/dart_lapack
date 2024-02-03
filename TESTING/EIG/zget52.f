      SUBROUTINE ZGET52( LEFT, N, A, LDA, B, LDB, E, LDE, ALPHA, BETA, WORK, RWORK, RESULT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      LOGICAL            LEFT
      int                LDA, LDB, LDE, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RESULT( 2 ), RWORK( * )
      COMPLEX*16         A( LDA, * ), ALPHA( * ), B( LDB, * ), BETA( * ), E( LDE, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      String             NORMAB, TRANS;
      int                J, JVEC
      DOUBLE PRECISION   ABMAX, ALFMAX, ANORM, BETMAX, BNORM, ENORM, ENRMER, ERRNRM, SAFMAX, SAFMIN, SCALE, TEMP1, ULP
      COMPLEX*16         ACOEFF, ALPHAI, BCOEFF, BETAI, X
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, ZLANGE
      EXTERNAL           DLAMCH, ZLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZGEMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCONJG, DIMAG, MAX
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   ABS1
*     ..
*     .. Statement Function definitions ..
      ABS1( X ) = ABS( DBLE( X ) ) + ABS( DIMAG( X ) )
*     ..
*     .. Executable Statements ..
*
      RESULT( 1 ) = ZERO
      RESULT( 2 ) = ZERO
      IF( N.LE.0 ) RETURN
*
      SAFMIN = DLAMCH( 'Safe minimum' )
      SAFMAX = ONE / SAFMIN
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' )
*
      IF( LEFT ) THEN
         TRANS = 'C'
         NORMAB = 'I'
      ELSE
         TRANS = 'N'
         NORMAB = 'O'
      END IF
*
*     Norm of A, B, and E:
*
      ANORM = MAX( ZLANGE( NORMAB, N, N, A, LDA, RWORK ), SAFMIN )
      BNORM = MAX( ZLANGE( NORMAB, N, N, B, LDB, RWORK ), SAFMIN )
      ENORM = MAX( ZLANGE( 'O', N, N, E, LDE, RWORK ), ULP )
      ALFMAX = SAFMAX / MAX( ONE, BNORM )
      BETMAX = SAFMAX / MAX( ONE, ANORM )
*
*     Compute error matrix.
*     Column i = ( b(i) A - a(i) B ) E(i) / max( |a(i) B|, |b(i) A| )
*
      DO 10 JVEC = 1, N
         ALPHAI = ALPHA( JVEC )
         BETAI = BETA( JVEC )
         ABMAX = MAX( ABS1( ALPHAI ), ABS1( BETAI ) )
         IF( ABS1( ALPHAI ).GT.ALFMAX .OR. ABS1( BETAI ).GT.BETMAX .OR. ABMAX.LT.ONE ) THEN
            SCALE = ONE / MAX( ABMAX, SAFMIN )
            ALPHAI = SCALE*ALPHAI
            BETAI = SCALE*BETAI
         END IF
         SCALE = ONE / MAX( ABS1( ALPHAI )*BNORM, ABS1( BETAI )*ANORM, SAFMIN )
         ACOEFF = SCALE*BETAI
         BCOEFF = SCALE*ALPHAI
         IF( LEFT ) THEN
            ACOEFF = DCONJG( ACOEFF )
            BCOEFF = DCONJG( BCOEFF )
         END IF
         CALL ZGEMV( TRANS, N, N, ACOEFF, A, LDA, E( 1, JVEC ), 1, CZERO, WORK( N*( JVEC-1 )+1 ), 1 )          CALL ZGEMV( TRANS, N, N, -BCOEFF, B, LDA, E( 1, JVEC ), 1, CONE, WORK( N*( JVEC-1 )+1 ), 1 )
   10 CONTINUE
*
      ERRNRM = ZLANGE( 'One', N, N, WORK, N, RWORK ) / ENORM
*
*     Compute RESULT(1)
*
      RESULT( 1 ) = ERRNRM / ULP
*
*     Normalization of E:
*
      ENRMER = ZERO
      DO 30 JVEC = 1, N
         TEMP1 = ZERO
         DO 20 J = 1, N
            TEMP1 = MAX( TEMP1, ABS1( E( J, JVEC ) ) )
   20    CONTINUE
         ENRMER = MAX( ENRMER, ABS( TEMP1-ONE ) )
   30 CONTINUE
*
*     Compute RESULT(2) : the normalization error in E.
*
      RESULT( 2 ) = ENRMER / ( DBLE( N )*ULP )
*
      RETURN
*
*     End of ZGET52
*
      END
