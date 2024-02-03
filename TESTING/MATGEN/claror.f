      SUBROUTINE CLAROR( SIDE, INIT, M, N, A, LDA, ISEED, X, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             INIT, SIDE;
      int                INFO, LDA, M, N;
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 );
      COMPLEX            A( LDA, * ), X( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, TOOSML
      const              ZERO = 0.0E+0, ONE = 1.0E+0, TOOSML = 1.0E-20 ;
      COMPLEX            CZERO, CONE
      const              CZERO = ( 0.0E+0, 0.0E+0 ), CONE = ( 1.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      int                IROW, ITYPE, IXFRM, J, JCOL, KBEG, NXFRM;
      REAL               FACTOR, XABS, XNORM
      COMPLEX            CSIGN, XNORMS
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SCNRM2
      COMPLEX            CLARND
      // EXTERNAL LSAME, SCNRM2, CLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMV, CGERC, CLACGV, CLASET, CSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CMPLX, CONJG
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
         CALL XERBLA( 'CLAROR', -INFO )
         RETURN
      END IF

      IF( ITYPE.EQ.1 ) THEN
         NXFRM = M
      } else {
         NXFRM = N
      END IF

      // Initialize A to the identity matrix if desired

      IF( LSAME( INIT, 'I' ) ) CALL CLASET( 'Full', M, N, CZERO, CONE, A, LDA )

      // If no rotation possible, still multiply by
      // a random complex number from the circle |x| = 1

       // 2)      Compute Rotation by computing Householder
               // Transformations H(2), H(3), ..., H(n).  Note that the
               // order in which they are computed is irrelevant.

      DO 40 J = 1, NXFRM
         X( J ) = CZERO
   40 CONTINUE

      DO 60 IXFRM = 2, NXFRM
         KBEG = NXFRM - IXFRM + 1

         // Generate independent normal( 0, 1 ) random numbers

         DO 50 J = KBEG, NXFRM
            X( J ) = CLARND( 3, ISEED )
   50    CONTINUE

         // Generate a Householder transformation from the random vector X

         XNORM = SCNRM2( IXFRM, X( KBEG ), 1 )
         XABS = ABS( X( KBEG ) )
         IF( XABS.NE.CZERO ) THEN
            CSIGN = X( KBEG ) / XABS
         } else {
            CSIGN = CONE
         END IF
         XNORMS = CSIGN*XNORM
         X( NXFRM+KBEG ) = -CSIGN
         FACTOR = XNORM*( XNORM+XABS )
         IF( ABS( FACTOR ).LT.TOOSML ) THEN
            INFO = 1
            CALL XERBLA( 'CLAROR', -INFO )
            RETURN
         } else {
            FACTOR = ONE / FACTOR
         END IF
         X( KBEG ) = X( KBEG ) + XNORMS

         // Apply Householder transformation to A

         IF( ITYPE.EQ.1 .OR. ITYPE.EQ.3 .OR. ITYPE.EQ.4 ) THEN

            // Apply H(k) on the left of A

            CALL CGEMV( 'C', IXFRM, N, CONE, A( KBEG, 1 ), LDA, X( KBEG ), 1, CZERO, X( 2*NXFRM+1 ), 1 )             CALL CGERC( IXFRM, N, -CMPLX( FACTOR ), X( KBEG ), 1, X( 2*NXFRM+1 ), 1, A( KBEG, 1 ), LDA )

         END IF

         IF( ITYPE.GE.2 .AND. ITYPE.LE.4 ) THEN

            // Apply H(k)* (or H(k)') on the right of A

            IF( ITYPE.EQ.4 ) THEN
               CALL CLACGV( IXFRM, X( KBEG ), 1 )
            END IF

            CALL CGEMV( 'N', M, IXFRM, CONE, A( 1, KBEG ), LDA, X( KBEG ), 1, CZERO, X( 2*NXFRM+1 ), 1 )             CALL CGERC( M, IXFRM, -CMPLX( FACTOR ), X( 2*NXFRM+1 ), 1, X( KBEG ), 1, A( 1, KBEG ), LDA )

         END IF
   60 CONTINUE

      X( 1 ) = CLARND( 3, ISEED )
      XABS = ABS( X( 1 ) )
      IF( XABS.NE.ZERO ) THEN
         CSIGN = X( 1 ) / XABS
      } else {
         CSIGN = CONE
      END IF
      X( 2*NXFRM ) = CSIGN

      // Scale the matrix A by D.

      IF( ITYPE.EQ.1 .OR. ITYPE.EQ.3 .OR. ITYPE.EQ.4 ) THEN
         DO 70 IROW = 1, M
            CALL CSCAL( N, CONJG( X( NXFRM+IROW ) ), A( IROW, 1 ), LDA )
   70    CONTINUE
      END IF

      IF( ITYPE.EQ.2 .OR. ITYPE.EQ.3 ) THEN
         DO 80 JCOL = 1, N
            CALL CSCAL( M, X( NXFRM+JCOL ), A( 1, JCOL ), 1 )
   80    CONTINUE
      END IF

      IF( ITYPE.EQ.4 ) THEN
         DO 90 JCOL = 1, N
            CALL CSCAL( M, CONJG( X( NXFRM+JCOL ) ), A( 1, JCOL ), 1 )
   90    CONTINUE
      END IF
      RETURN

      // End of CLAROR

      }
