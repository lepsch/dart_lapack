      SUBROUTINE CGGHRD( COMPQ, COMPZ, N, ILO, IHI, A, LDA, B, LDB, Q, LDQ, Z, LDZ, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             COMPQ, COMPZ;
      int                IHI, ILO, INFO, LDA, LDB, LDQ, LDZ, N;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * ), Q( LDQ, * ), Z( LDZ, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            CONE, CZERO
      const              CONE = ( 1.0E+0, 0.0E+0 ), CZERO = ( 0.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               ILQ, ILZ;
      int                ICOMPQ, ICOMPZ, JCOL, JROW;
      REAL               C
      COMPLEX            CTEMP, S
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLARTG, CLASET, CROT, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG, MAX
      // ..
      // .. Executable Statements ..

      // Decode COMPQ

      IF( LSAME( COMPQ, 'N' ) ) THEN
         ILQ = .FALSE.
         ICOMPQ = 1
      ELSE IF( LSAME( COMPQ, 'V' ) ) THEN
         ILQ = .TRUE.
         ICOMPQ = 2
      ELSE IF( LSAME( COMPQ, 'I' ) ) THEN
         ILQ = .TRUE.
         ICOMPQ = 3
      } else {
         ICOMPQ = 0
      END IF

      // Decode COMPZ

      IF( LSAME( COMPZ, 'N' ) ) THEN
         ILZ = .FALSE.
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ILZ = .TRUE.
         ICOMPZ = 2
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ILZ = .TRUE.
         ICOMPZ = 3
      } else {
         ICOMPZ = 0
      END IF

      // Test the input parameters.

      INFO = 0
      IF( ICOMPQ.LE.0 ) THEN
         INFO = -1
      ELSE IF( ICOMPZ.LE.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( ILO.LT.1 ) THEN
         INFO = -4
      ELSE IF( IHI.GT.N .OR. IHI.LT.ILO-1 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( ( ILQ .AND. LDQ.LT.N ) .OR. LDQ.LT.1 ) THEN
         INFO = -11
      ELSE IF( ( ILZ .AND. LDZ.LT.N ) .OR. LDZ.LT.1 ) THEN
         INFO = -13
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGGHRD', -INFO )
         RETURN
      END IF

      // Initialize Q and Z if desired.

      IF( ICOMPQ.EQ.3 ) CALL CLASET( 'Full', N, N, CZERO, CONE, Q, LDQ )       IF( ICOMPZ.EQ.3 ) CALL CLASET( 'Full', N, N, CZERO, CONE, Z, LDZ )

      // Quick return if possible

      IF( N.LE.1 ) RETURN

      // Zero out lower triangle of B

      DO 20 JCOL = 1, N - 1
         DO 10 JROW = JCOL + 1, N
            B( JROW, JCOL ) = CZERO
   10    CONTINUE
   20 CONTINUE

      // Reduce A and B

      DO 40 JCOL = ILO, IHI - 2

         DO 30 JROW = IHI, JCOL + 2, -1

            // Step 1: rotate rows JROW-1, JROW to kill A(JROW,JCOL)

            CTEMP = A( JROW-1, JCOL )
            CALL CLARTG( CTEMP, A( JROW, JCOL ), C, S, A( JROW-1, JCOL ) )
            A( JROW, JCOL ) = CZERO
            CALL CROT( N-JCOL, A( JROW-1, JCOL+1 ), LDA, A( JROW, JCOL+1 ), LDA, C, S )             CALL CROT( N+2-JROW, B( JROW-1, JROW-1 ), LDB, B( JROW, JROW-1 ), LDB, C, S )             IF( ILQ ) CALL CROT( N, Q( 1, JROW-1 ), 1, Q( 1, JROW ), 1, C, CONJG( S ) )

            // Step 2: rotate columns JROW, JROW-1 to kill B(JROW,JROW-1)

            CTEMP = B( JROW, JROW )
            CALL CLARTG( CTEMP, B( JROW, JROW-1 ), C, S, B( JROW, JROW ) )
            B( JROW, JROW-1 ) = CZERO
            CALL CROT( IHI, A( 1, JROW ), 1, A( 1, JROW-1 ), 1, C, S )
            CALL CROT( JROW-1, B( 1, JROW ), 1, B( 1, JROW-1 ), 1, C, S )             IF( ILZ ) CALL CROT( N, Z( 1, JROW ), 1, Z( 1, JROW-1 ), 1, C, S )
   30    CONTINUE
   40 CONTINUE

      RETURN

      // End of CGGHRD

      }
