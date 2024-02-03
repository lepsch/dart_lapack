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
      IF( N.EQ.0 .OR. M.EQ.0 ) RETURN

      ITYPE = 0
      if ( LSAME( SIDE, 'L' ) ) {
         ITYPE = 1
      } else if ( LSAME( SIDE, 'R' ) ) {
         ITYPE = 2
      } else if ( LSAME( SIDE, 'C' ) .OR. LSAME( SIDE, 'T' ) ) {
         ITYPE = 3
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
         xerbla('DLAROR', -INFO );
         RETURN
      }

      if ( ITYPE.EQ.1 ) {
         NXFRM = M
      } else {
         NXFRM = N
      }

      // Initialize A to the identity matrix if desired

      IF( LSAME( INIT, 'I' ) ) CALL DLASET( 'Full', M, N, ZERO, ONE, A, LDA )

      // If no rotation possible, multiply by random +/-1

      // Compute rotation by computing Householder transformations
      // H(2), H(3), ..., H(nhouse)

      DO 10 J = 1, NXFRM
         X( J ) = ZERO
   10 CONTINUE

      DO 30 IXFRM = 2, NXFRM
         KBEG = NXFRM - IXFRM + 1

         // Generate independent normal( 0, 1 ) random numbers

         DO 20 J = KBEG, NXFRM
            X( J ) = DLARND( 3, ISEED )
   20    CONTINUE

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

         if ( ITYPE.EQ.1 .OR. ITYPE.EQ.3 ) {

            // Apply H(k) from the left.

            dgemv('T', IXFRM, N, ONE, A( KBEG, 1 ), LDA, X( KBEG ), 1, ZERO, X( 2*NXFRM+1 ), 1 )             CALL DGER( IXFRM, N, -FACTOR, X( KBEG ), 1, X( 2*NXFRM+1 ), 1, A( KBEG, 1 ), LDA );

         }

         if ( ITYPE.EQ.2 .OR. ITYPE.EQ.3 ) {

            // Apply H(k) from the right.

            dgemv('N', M, IXFRM, ONE, A( 1, KBEG ), LDA, X( KBEG ), 1, ZERO, X( 2*NXFRM+1 ), 1 )             CALL DGER( M, IXFRM, -FACTOR, X( 2*NXFRM+1 ), 1, X( KBEG ), 1, A( 1, KBEG ), LDA );

         }
   30 CONTINUE

      X( 2*NXFRM ) = SIGN( ONE, DLARND( 3, ISEED ) )

      // Scale the matrix A by D.

      if ( ITYPE.EQ.1 .OR. ITYPE.EQ.3 ) {
         DO 40 IROW = 1, M
            dscal(N, X( NXFRM+IROW ), A( IROW, 1 ), LDA );
   40    CONTINUE
      }

      if ( ITYPE.EQ.2 .OR. ITYPE.EQ.3 ) {
         DO 50 JCOL = 1, N
            dscal(M, X( NXFRM+JCOL ), A( 1, JCOL ), 1 );
   50    CONTINUE
      }
      RETURN

      // End of DLAROR

      }
