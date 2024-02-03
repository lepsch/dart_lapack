      SUBROUTINE ZUNT01( ROWCOL, M, N, U, LDU, WORK, LWORK, RWORK, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             ROWCOL;
      int                LDU, LWORK, M, N
      DOUBLE PRECISION   RESID
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         U( LDU, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      String             TRANSU;
      int                I, J, K, LDWORK, MNMIN
      DOUBLE PRECISION   EPS
      COMPLEX*16         TMP, ZDUM
*     ..
*     .. External Functions ..
      bool               LSAME;
      DOUBLE PRECISION   DLAMCH, ZLANSY
      COMPLEX*16         ZDOTC
      EXTERNAL           LSAME, DLAMCH, ZLANSY, ZDOTC
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZHERK, ZLASET
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, DIMAG, MAX, MIN
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )
*     ..
*     .. Executable Statements ..
*
      RESID = ZERO
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 ) RETURN
*
      EPS = DLAMCH( 'Precision' )
      IF( M.LT.N .OR. ( M.EQ.N .AND. LSAME( ROWCOL, 'R' ) ) ) THEN
         TRANSU = 'N'
         K = N
      ELSE
         TRANSU = 'C'
         K = M
      END IF
      MNMIN = MIN( M, N )
*
      IF( ( MNMIN+1 )*MNMIN.LE.LWORK ) THEN
         LDWORK = MNMIN
      ELSE
         LDWORK = 0
      END IF
      IF( LDWORK.GT.0 ) THEN
*
*        Compute I - U*U' or I - U'*U.
*
         CALL ZLASET( 'Upper', MNMIN, MNMIN, DCMPLX( ZERO ), DCMPLX( ONE ), WORK, LDWORK )          CALL ZHERK( 'Upper', TRANSU, MNMIN, K, -ONE, U, LDU, ONE, WORK, LDWORK )
*
*        Compute norm( I - U*U' ) / ( K * EPS ) .
*
         RESID = ZLANSY( '1', 'Upper', MNMIN, WORK, LDWORK, RWORK )
         RESID = ( RESID / DBLE( K ) ) / EPS
      ELSE IF( TRANSU.EQ.'C' ) THEN
*
*        Find the maximum element in abs( I - U'*U ) / ( m * EPS )
*
         DO 20 J = 1, N
            DO 10 I = 1, J
               IF( I.NE.J ) THEN
                  TMP = ZERO
               ELSE
                  TMP = ONE
               END IF
               TMP = TMP - ZDOTC( M, U( 1, I ), 1, U( 1, J ), 1 )
               RESID = MAX( RESID, CABS1( TMP ) )
   10       CONTINUE
   20    CONTINUE
         RESID = ( RESID / DBLE( M ) ) / EPS
      ELSE
*
*        Find the maximum element in abs( I - U*U' ) / ( n * EPS )
*
         DO 40 J = 1, M
            DO 30 I = 1, J
               IF( I.NE.J ) THEN
                  TMP = ZERO
               ELSE
                  TMP = ONE
               END IF
               TMP = TMP - ZDOTC( N, U( J, 1 ), LDU, U( I, 1 ), LDU )
               RESID = MAX( RESID, CABS1( TMP ) )
   30       CONTINUE
   40    CONTINUE
         RESID = ( RESID / DBLE( N ) ) / EPS
      END IF
      RETURN
*
*     End of ZUNT01
*
      END
