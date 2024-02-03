      SUBROUTINE CUNT01( ROWCOL, M, N, U, LDU, WORK, LWORK, RWORK, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          ROWCOL
      INTEGER            LDU, LWORK, M, N
      REAL               RESID
*     ..
*     .. Array Arguments ..
      REAL               RWORK( * )
      COMPLEX            U( LDU, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      CHARACTER          TRANSU
      INTEGER            I, J, K, LDWORK, MNMIN
      REAL               EPS
      COMPLEX            TMP, ZDUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               CLANSY, SLAMCH
      COMPLEX            CDOTC
      EXTERNAL           LSAME, CLANSY, SLAMCH, CDOTC
*     ..
*     .. External Subroutines ..
      EXTERNAL           CHERK, CLASET
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, AIMAG, CMPLX, MAX, MIN, REAL
*     ..
*     .. Statement Functions ..
      REAL               CABS1
*     ..
*     .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )
*     ..
*     .. Executable Statements ..
*
      RESID = ZERO
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 ) RETURN
*
      EPS = SLAMCH( 'Precision' )
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
         CALL CLASET( 'Upper', MNMIN, MNMIN, CMPLX( ZERO ), CMPLX( ONE ), WORK, LDWORK )          CALL CHERK( 'Upper', TRANSU, MNMIN, K, -ONE, U, LDU, ONE, WORK, LDWORK )
*
*        Compute norm( I - U*U' ) / ( K * EPS ) .
*
         RESID = CLANSY( '1', 'Upper', MNMIN, WORK, LDWORK, RWORK )
         RESID = ( RESID / REAL( K ) ) / EPS
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
               TMP = TMP - CDOTC( M, U( 1, I ), 1, U( 1, J ), 1 )
               RESID = MAX( RESID, CABS1( TMP ) )
   10       CONTINUE
   20    CONTINUE
         RESID = ( RESID / REAL( M ) ) / EPS
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
               TMP = TMP - CDOTC( N, U( J, 1 ), LDU, U( I, 1 ), LDU )
               RESID = MAX( RESID, CABS1( TMP ) )
   30       CONTINUE
   40    CONTINUE
         RESID = ( RESID / REAL( N ) ) / EPS
      END IF
      RETURN
*
*     End of CUNT01
*
      END
