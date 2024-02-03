      SUBROUTINE ZLAROR( SIDE, INIT, M, N, A, LDA, ISEED, X, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             INIT, SIDE;
      int                INFO, LDA, M, N;
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 );
      COMPLEX*16         A( LDA, * ), X( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TOOSML;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TOOSML = 1.0D-20 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) )
      // ..
      // .. Local Scalars ..
      int                IROW, ITYPE, IXFRM, J, JCOL, KBEG, NXFRM;
      double             FACTOR, XABS, XNORM;
      COMPLEX*16         CSIGN, XNORMS
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DZNRM2;
      COMPLEX*16         ZLARND
      // EXTERNAL LSAME, DZNRM2, ZLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZGEMV, ZGERC, ZLACGV, ZLASET, ZSCAL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DCMPLX, DCONJG
      // ..
      // .. Executable Statements ..

      INFO = 0
      IF( N.EQ.0 .OR. M.EQ.0 ) RETURN

      ITYPE = 0
      IF( LSAME( SIDE, 'L' ) ) THEN
         ITYPE = 1
      ELSE IF( LSAME( SIDE, 'R' ) ) THEN
         ITYPE = 2
      ELSE IF( LSAME( SIDE, 'C' ) ) THEN
         ITYPE = 3
      ELSE IF( LSAME( SIDE, 'T' ) ) THEN
         ITYPE = 4
      END IF

      // Check for argument errors.

      IF( ITYPE.EQ.0 ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 .OR. ( ITYPE.EQ.3 .AND. N.NE.M ) ) THEN
         INFO = -4
      ELSE IF( LDA.LT.M ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZLAROR', -INFO )
         RETURN
      END IF

      IF( ITYPE.EQ.1 ) THEN
         NXFRM = M
      ELSE
         NXFRM = N
      END IF

      // Initialize A to the identity matrix if desired

      IF( LSAME( INIT, 'I' ) ) CALL ZLASET( 'Full', M, N, CZERO, CONE, A, LDA )

      // If no rotation possible, still multiply by
      // a random complex number from the circle |x| = 1

       // 2)      Compute Rotation by computing Householder
               // Transformations H(2), H(3), ..., H(n).  Note that the
               // order in which they are computed is irrelevant.

      DO 10 J = 1, NXFRM
         X( J ) = CZERO
   10 CONTINUE

      DO 30 IXFRM = 2, NXFRM
         KBEG = NXFRM - IXFRM + 1

         // Generate independent normal( 0, 1 ) random numbers

         DO 20 J = KBEG, NXFRM
            X( J ) = ZLARND( 3, ISEED )
   20    CONTINUE

         // Generate a Householder transformation from the random vector X

         XNORM = DZNRM2( IXFRM, X( KBEG ), 1 )
         XABS = ABS( X( KBEG ) )
         IF( XABS.NE.CZERO ) THEN
            CSIGN = X( KBEG ) / XABS
         ELSE
            CSIGN = CONE
         END IF
         XNORMS = CSIGN*XNORM
         X( NXFRM+KBEG ) = -CSIGN
         FACTOR = XNORM*( XNORM+XABS )
         IF( ABS( FACTOR ).LT.TOOSML ) THEN
            INFO = 1
            CALL XERBLA( 'ZLAROR', -INFO )
            RETURN
         ELSE
            FACTOR = ONE / FACTOR
         END IF
         X( KBEG ) = X( KBEG ) + XNORMS

         // Apply Householder transformation to A

         IF( ITYPE.EQ.1 .OR. ITYPE.EQ.3 .OR. ITYPE.EQ.4 ) THEN

            // Apply H(k) on the left of A

            CALL ZGEMV( 'C', IXFRM, N, CONE, A( KBEG, 1 ), LDA, X( KBEG ), 1, CZERO, X( 2*NXFRM+1 ), 1 )             CALL ZGERC( IXFRM, N, -DCMPLX( FACTOR ), X( KBEG ), 1, X( 2*NXFRM+1 ), 1, A( KBEG, 1 ), LDA )

         END IF

         IF( ITYPE.GE.2 .AND. ITYPE.LE.4 ) THEN

            // Apply H(k)* (or H(k)') on the right of A

            IF( ITYPE.EQ.4 ) THEN
               CALL ZLACGV( IXFRM, X( KBEG ), 1 )
            END IF

            CALL ZGEMV( 'N', M, IXFRM, CONE, A( 1, KBEG ), LDA, X( KBEG ), 1, CZERO, X( 2*NXFRM+1 ), 1 )             CALL ZGERC( M, IXFRM, -DCMPLX( FACTOR ), X( 2*NXFRM+1 ), 1, X( KBEG ), 1, A( 1, KBEG ), LDA )

         END IF
   30 CONTINUE

      X( 1 ) = ZLARND( 3, ISEED )
      XABS = ABS( X( 1 ) )
      IF( XABS.NE.ZERO ) THEN
         CSIGN = X( 1 ) / XABS
      ELSE
         CSIGN = CONE
      END IF
      X( 2*NXFRM ) = CSIGN

      // Scale the matrix A by D.

      IF( ITYPE.EQ.1 .OR. ITYPE.EQ.3 .OR. ITYPE.EQ.4 ) THEN
         DO 40 IROW = 1, M
            CALL ZSCAL( N, DCONJG( X( NXFRM+IROW ) ), A( IROW, 1 ), LDA )
   40    CONTINUE
      END IF

      IF( ITYPE.EQ.2 .OR. ITYPE.EQ.3 ) THEN
         DO 50 JCOL = 1, N
            CALL ZSCAL( M, X( NXFRM+JCOL ), A( 1, JCOL ), 1 )
   50    CONTINUE
      END IF

      IF( ITYPE.EQ.4 ) THEN
         DO 60 JCOL = 1, N
            CALL ZSCAL( M, DCONJG( X( NXFRM+JCOL ) ), A( 1, JCOL ), 1 )
   60    CONTINUE
      END IF
      RETURN

      // End of ZLAROR

      END
