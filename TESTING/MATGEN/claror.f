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
      if (N == 0 .OR. M == 0) RETURN;

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

      if ( ITYPE == 0 ) {
         INFO = -1
      } else if ( M.LT.0 ) {
         INFO = -3
      } else if ( N.LT.0 .OR. ( ITYPE == 3 .AND. N.NE.M ) ) {
         INFO = -4
      } else if ( LDA.LT.M ) {
         INFO = -6
      }
      if ( INFO.NE.0 ) {
         xerbla('CLAROR', -INFO );
         RETURN
      }

      if ( ITYPE == 1 ) {
         NXFRM = M
      } else {
         NXFRM = N
      }

      // Initialize A to the identity matrix if desired

      IF( LSAME( INIT, 'I' ) ) CALL CLASET( 'Full', M, N, CZERO, CONE, A, LDA )

      // If no rotation possible, still multiply by
      // a random complex number from the circle |x| = 1

       // 2)      Compute Rotation by computing Householder
               // Transformations H(2), H(3), ..., H(n).  Note that the
               // order in which they are computed is irrelevant.

      for (J = 1; J <= NXFRM; J++) { // 40
         X( J ) = CZERO
      } // 40

      for (IXFRM = 2; IXFRM <= NXFRM; IXFRM++) { // 60
         KBEG = NXFRM - IXFRM + 1

         // Generate independent normal( 0, 1 ) random numbers

         for (J = KBEG; J <= NXFRM; J++) { // 50
            X( J ) = CLARND( 3, ISEED )
         } // 50

         // Generate a Householder transformation from the random vector X

         XNORM = SCNRM2( IXFRM, X( KBEG ), 1 )
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
            xerbla('CLAROR', -INFO );
            RETURN
         } else {
            FACTOR = ONE / FACTOR
         }
         X( KBEG ) = X( KBEG ) + XNORMS

         // Apply Householder transformation to A

         if ( ITYPE == 1 .OR. ITYPE == 3 .OR. ITYPE == 4 ) {

            // Apply H(k) on the left of A

            cgemv('C', IXFRM, N, CONE, A( KBEG, 1 ), LDA, X( KBEG ), 1, CZERO, X( 2*NXFRM+1 ), 1 );
            cgerc(IXFRM, N, -CMPLX( FACTOR ), X( KBEG ), 1, X( 2*NXFRM+1 ), 1, A( KBEG, 1 ), LDA );

         }

         if ( ITYPE.GE.2 .AND. ITYPE.LE.4 ) {

            // Apply H(k)* (or H(k)') on the right of A

            if ( ITYPE == 4 ) {
               clacgv(IXFRM, X( KBEG ), 1 );
            }

            cgemv('N', M, IXFRM, CONE, A( 1, KBEG ), LDA, X( KBEG ), 1, CZERO, X( 2*NXFRM+1 ), 1 );
            cgerc(M, IXFRM, -CMPLX( FACTOR ), X( 2*NXFRM+1 ), 1, X( KBEG ), 1, A( 1, KBEG ), LDA );

         }
      } // 60

      X( 1 ) = CLARND( 3, ISEED )
      XABS = ABS( X( 1 ) )
      if ( XABS.NE.ZERO ) {
         CSIGN = X( 1 ) / XABS
      } else {
         CSIGN = CONE
      }
      X( 2*NXFRM ) = CSIGN

      // Scale the matrix A by D.

      if ( ITYPE == 1 .OR. ITYPE == 3 .OR. ITYPE == 4 ) {
         for (IROW = 1; IROW <= M; IROW++) { // 70
            cscal(N, CONJG( X( NXFRM+IROW ) ), A( IROW, 1 ), LDA );
         } // 70
      }

      if ( ITYPE == 2 .OR. ITYPE == 3 ) {
         for (JCOL = 1; JCOL <= N; JCOL++) { // 80
            cscal(M, X( NXFRM+JCOL ), A( 1, JCOL ), 1 );
         } // 80
      }

      if ( ITYPE == 4 ) {
         for (JCOL = 1; JCOL <= N; JCOL++) { // 90
            cscal(M, CONJG( X( NXFRM+JCOL ) ), A( 1, JCOL ), 1 );
         } // 90
      }
      RETURN

      // End of CLAROR

      }
