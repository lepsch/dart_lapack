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
      const              ZERO = 0.0D+0, ONE = 1.0D+0, TOOSML = 1.0D-20 ;
      COMPLEX*16         CZERO, CONE
      const              CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) ;
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
      if ( LSAME( SIDE, 'L' ) ) {
         ITYPE = 1
      } else if ( LSAME( SIDE, 'R' ) ) {
         ITYPE = 2
      } else if ( LSAME( SIDE, 'C' ) ) {
         ITYPE = 3
      } else if ( LSAME( SIDE, 'T' ) ) {
         ITYPE = 4
      }

      // Check for argument errors.

      if ( ITYPE.EQ.0 ) {
         INFO = -1
      } else if ( M.LT.0 ) {
         INFO = -3
      } else if ( N.LT.0 .OR. ( ITYPE.EQ.3 .AND. N.NE.M ) ) {
         INFO = -4
      } else if ( LDA.LT.M ) {
         INFO = -6
      }
      if ( INFO.NE.0 ) {
         xerbla('ZLAROR', -INFO );
         RETURN
      }

      if ( ITYPE.EQ.1 ) {
         NXFRM = M
      } else {
         NXFRM = N
      }

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
         if ( XABS.NE.CZERO ) {
            CSIGN = X( KBEG ) / XABS
         } else {
            CSIGN = CONE
         }
         XNORMS = CSIGN*XNORM
         X( NXFRM+KBEG ) = -CSIGN
         FACTOR = XNORM*( XNORM+XABS )
         if ( ABS( FACTOR ).LT.TOOSML ) {
            INFO = 1
            xerbla('ZLAROR', -INFO );
            RETURN
         } else {
            FACTOR = ONE / FACTOR
         }
         X( KBEG ) = X( KBEG ) + XNORMS

         // Apply Householder transformation to A

         if ( ITYPE.EQ.1 .OR. ITYPE.EQ.3 .OR. ITYPE.EQ.4 ) {

            // Apply H(k) on the left of A

            zgemv('C', IXFRM, N, CONE, A( KBEG, 1 ), LDA, X( KBEG ), 1, CZERO, X( 2*NXFRM+1 ), 1 )             CALL ZGERC( IXFRM, N, -DCMPLX( FACTOR ), X( KBEG ), 1, X( 2*NXFRM+1 ), 1, A( KBEG, 1 ), LDA );

         }

         if ( ITYPE.GE.2 .AND. ITYPE.LE.4 ) {

            // Apply H(k)* (or H(k)') on the right of A

            if ( ITYPE.EQ.4 ) {
               zlacgv(IXFRM, X( KBEG ), 1 );
            }

            zgemv('N', M, IXFRM, CONE, A( 1, KBEG ), LDA, X( KBEG ), 1, CZERO, X( 2*NXFRM+1 ), 1 )             CALL ZGERC( M, IXFRM, -DCMPLX( FACTOR ), X( 2*NXFRM+1 ), 1, X( KBEG ), 1, A( 1, KBEG ), LDA );

         }
   30 CONTINUE

      X( 1 ) = ZLARND( 3, ISEED )
      XABS = ABS( X( 1 ) )
      if ( XABS.NE.ZERO ) {
         CSIGN = X( 1 ) / XABS
      } else {
         CSIGN = CONE
      }
      X( 2*NXFRM ) = CSIGN

      // Scale the matrix A by D.

      if ( ITYPE.EQ.1 .OR. ITYPE.EQ.3 .OR. ITYPE.EQ.4 ) {
         DO 40 IROW = 1, M
            zscal(N, DCONJG( X( NXFRM+IROW ) ), A( IROW, 1 ), LDA );
   40    CONTINUE
      }

      if ( ITYPE.EQ.2 .OR. ITYPE.EQ.3 ) {
         DO 50 JCOL = 1, N
            zscal(M, X( NXFRM+JCOL ), A( 1, JCOL ), 1 );
   50    CONTINUE
      }

      if ( ITYPE.EQ.4 ) {
         DO 60 JCOL = 1, N
            zscal(M, DCONJG( X( NXFRM+JCOL ) ), A( 1, JCOL ), 1 );
   60    CONTINUE
      }
      RETURN

      // End of ZLAROR

      }
