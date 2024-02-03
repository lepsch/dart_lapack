      SUBROUTINE ZPTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             COMPZ;
      int                INFO, LDZ, N;
      // ..
      // .. Array Arguments ..
      double             D( * ), E( * ), WORK( * );
      COMPLEX*16         Z( LDZ, * )
      // ..

*  ====================================================================

      // .. Parameters ..
      COMPLEX*16         CZERO, CONE
      const              CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) ;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL DPTTRF, XERBLA, ZBDSQR, ZLASET
      // ..
      // .. Local Arrays ..
      COMPLEX*16         C( 1, 1 ), VT( 1, 1 )
      // ..
      // .. Local Scalars ..
      int                I, ICOMPZ, NRU;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0

      IF( LSAME( COMPZ, 'N' ) ) THEN
         ICOMPZ = 0
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ICOMPZ = 2
      } else {
         ICOMPZ = -1
      END IF
      IF( ICOMPZ.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1, N ) ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZPTEQR', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      IF( N.EQ.1 ) THEN
         IF( ICOMPZ.GT.0 ) Z( 1, 1 ) = CONE
         RETURN
      END IF
      IF( ICOMPZ.EQ.2 ) CALL ZLASET( 'Full', N, N, CZERO, CONE, Z, LDZ )

      // Call DPTTRF to factor the matrix.

      CALL DPTTRF( N, D, E, INFO )
      IF( INFO.NE.0 ) RETURN
      DO 10 I = 1, N
         D( I ) = SQRT( D( I ) )
   10 CONTINUE
      DO 20 I = 1, N - 1
         E( I ) = E( I )*D( I )
   20 CONTINUE

      // Call ZBDSQR to compute the singular values/vectors of the
      // bidiagonal factor.

      IF( ICOMPZ.GT.0 ) THEN
         NRU = N
      } else {
         NRU = 0
      END IF
      CALL ZBDSQR( 'Lower', N, 0, NRU, 0, D, E, VT, 1, Z, LDZ, C, 1, WORK, INFO )

      // Square the singular values.

      IF( INFO.EQ.0 ) THEN
         DO 30 I = 1, N
            D( I ) = D( I )*D( I )
   30    CONTINUE
      } else {
         INFO = N + INFO
      END IF

      RETURN

      // End of ZPTEQR

      }
