      SUBROUTINE DLAROR( SIDE, INIT, M, N, A, LDA, ISEED, X, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             INIT, SIDE;
      int                INFO, LDA, M, N;
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 );
      double             A( LDA, * ), X( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TOOSML;
      const              ZERO = 0.0D+0, ONE = 1.0D+0, TOOSML = 1.0D-20 ;
      // ..
      // .. Local Scalars ..
      int                IROW, ITYPE, IXFRM, J, JCOL, KBEG, NXFRM;
      double             FACTOR, XNORM, XNORMS;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLARND, DNRM2;
      // EXTERNAL LSAME, DLARND, DNRM2
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMV, DGER, DLASET, DSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, SIGN
      // ..
      // .. Executable Statements ..

      INFO = 0
      if (N == 0 .OR. M == 0) RETURN;

      ITYPE = 0
      if ( LSAME( SIDE, 'L' ) ) {
         ITYPE = 1
      } else if ( LSAME( SIDE, 'R' ) ) {
         ITYPE = 2
      } else if ( LSAME( SIDE, 'C' ) .OR. LSAME( SIDE, 'T' ) ) {
         ITYPE = 3
      }

      // Check for argument errors.

      if ( ITYPE == 0 ) {
         INFO = -1
      } else if ( M.LT.0 ) {
         INFO = -3
      } else if ( N.LT.0 .OR. ( ITYPE == 3 && N != M ) ) {
         INFO = -4
      } else if ( LDA.LT.M ) {
         INFO = -6
      }
      if ( INFO != 0 ) {
         xerbla('DLAROR', -INFO );
         RETURN
      }

      if ( ITYPE == 1 ) {
         NXFRM = M
      } else {
         NXFRM = N
      }

      // Initialize A to the identity matrix if desired

      IF( LSAME( INIT, 'I' ) ) CALL DLASET( 'Full', M, N, ZERO, ONE, A, LDA )

      // If no rotation possible, multiply by random +/-1

      // Compute rotation by computing Householder transformations
      // H(2), H(3), ..., H(nhouse)

      for (J = 1; J <= NXFRM; J++) { // 10
         X( J ) = ZERO
      } // 10

      for (IXFRM = 2; IXFRM <= NXFRM; IXFRM++) { // 30
         KBEG = NXFRM - IXFRM + 1

         // Generate independent normal( 0, 1 ) random numbers

         for (J = KBEG; J <= NXFRM; J++) { // 20
            X( J ) = DLARND( 3, ISEED )
         } // 20

         // Generate a Householder transformation from the random vector X

         XNORM = DNRM2( IXFRM, X( KBEG ), 1 )
         XNORMS = SIGN( XNORM, X( KBEG ) )
         X( KBEG+NXFRM ) = SIGN( ONE, -X( KBEG ) )
         FACTOR = XNORMS*( XNORMS+X( KBEG ) )
         if ( ABS( FACTOR ).LT.TOOSML ) {
            INFO = 1
            xerbla('DLAROR', INFO );
            RETURN
         } else {
            FACTOR = ONE / FACTOR
         }
         X( KBEG ) = X( KBEG ) + XNORMS

         // Apply Householder transformation to A

         if ( ITYPE == 1 .OR. ITYPE == 3 ) {

            // Apply H(k) from the left.

            dgemv('T', IXFRM, N, ONE, A( KBEG, 1 ), LDA, X( KBEG ), 1, ZERO, X( 2*NXFRM+1 ), 1 );
            dger(IXFRM, N, -FACTOR, X( KBEG ), 1, X( 2*NXFRM+1 ), 1, A( KBEG, 1 ), LDA );

         }

         if ( ITYPE == 2 .OR. ITYPE == 3 ) {

            // Apply H(k) from the right.

            dgemv('N', M, IXFRM, ONE, A( 1, KBEG ), LDA, X( KBEG ), 1, ZERO, X( 2*NXFRM+1 ), 1 );
            dger(M, IXFRM, -FACTOR, X( 2*NXFRM+1 ), 1, X( KBEG ), 1, A( 1, KBEG ), LDA );

         }
      } // 30

      X( 2*NXFRM ) = SIGN( ONE, DLARND( 3, ISEED ) )

      // Scale the matrix A by D.

      if ( ITYPE == 1 .OR. ITYPE == 3 ) {
         for (IROW = 1; IROW <= M; IROW++) { // 40
            dscal(N, X( NXFRM+IROW ), A( IROW, 1 ), LDA );
         } // 40
      }

      if ( ITYPE == 2 .OR. ITYPE == 3 ) {
         for (JCOL = 1; JCOL <= N; JCOL++) { // 50
            dscal(M, X( NXFRM+JCOL ), A( 1, JCOL ), 1 );
         } // 50
      }
      RETURN

      // End of DLAROR

      }
