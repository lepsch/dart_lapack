      void zlaror(SIDE, INIT, M, N, A, LDA, ISEED, X, INFO ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             INIT, SIDE;
      int                INFO, LDA, M, N;
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 );
      Complex         A( LDA, * ), X( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TOOSML;
      const              ZERO = 0.0, ONE = 1.0, TOOSML = 1.0e-20 ;
      Complex         CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      int                IROW, ITYPE, IXFRM, J, JCOL, KBEG, NXFRM;
      double             FACTOR, XABS, XNORM;
      Complex         CSIGN, XNORMS;
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      //- double             DZNRM2;
      //- Complex         ZLARND;
      // EXTERNAL LSAME, DZNRM2, ZLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZGEMV, ZGERC, ZLACGV, ZLASET, ZSCAL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DCMPLX, DCONJG
      // ..
      // .. Executable Statements ..

      INFO = 0;
      if (N == 0 || M == 0) return;

      ITYPE = 0;
      if ( LSAME( SIDE, 'L' ) ) {
         ITYPE = 1;
      } else if ( LSAME( SIDE, 'R' ) ) {
         ITYPE = 2;
      } else if ( LSAME( SIDE, 'C' ) ) {
         ITYPE = 3;
      } else if ( LSAME( SIDE, 'T' ) ) {
         ITYPE = 4;
      }

      // Check for argument errors.

      if ( ITYPE == 0 ) {
         INFO = -1;
      } else if ( M < 0 ) {
         INFO = -3;
      } else if ( N < 0 || ( ITYPE == 3 && N != M ) ) {
         INFO = -4;
      } else if ( LDA < M ) {
         INFO = -6;
      }
      if ( INFO != 0 ) {
         xerbla('ZLAROR', -INFO );
         return;
      }

      if ( ITYPE == 1 ) {
         NXFRM = M;
      } else {
         NXFRM = N;
      }

      // Initialize A to the identity matrix if desired

      if( LSAME( INIT, 'I' ) ) zlaset( 'Full', M, N, CZERO, CONE, A, LDA );

      // If no rotation possible, still multiply by
      // a random complex number from the circle |x| = 1

       // 2)      Compute Rotation by computing Householder
               // Transformations H(2), H(3), ..., H(n).  Note that the
               // order in which they are computed is irrelevant.

      for (J = 1; J <= NXFRM; J++) { // 10
         X( J ) = CZERO;
      } // 10

      for (IXFRM = 2; IXFRM <= NXFRM; IXFRM++) { // 30
         KBEG = NXFRM - IXFRM + 1;

         // Generate independent normal( 0, 1 ) random numbers

         for (J = KBEG; J <= NXFRM; J++) { // 20
            X( J ) = ZLARND( 3, ISEED );
         } // 20

         // Generate a Householder transformation from the random vector X

         XNORM = DZNRM2( IXFRM, X( KBEG ), 1 );
         XABS = ( X( KBEG ) ).abs();
         if ( XABS != CZERO ) {
            CSIGN = X( KBEG ) / XABS;
         } else {
            CSIGN = CONE;
         }
         XNORMS = CSIGN*XNORM;
         X( NXFRM+KBEG ) = -CSIGN;
         FACTOR = XNORM*( XNORM+XABS );
         if ( ( FACTOR ).abs() < TOOSML ) {
            INFO = 1;
            xerbla('ZLAROR', -INFO );
            return;
         } else {
            FACTOR = ONE / FACTOR;
         }
         X( KBEG ) = X( KBEG ) + XNORMS;

         // Apply Householder transformation to A

         if ( ITYPE == 1 || ITYPE == 3 || ITYPE == 4 ) {

            // Apply H(k) on the left of A

            zgemv('C', IXFRM, N, CONE, A( KBEG, 1 ), LDA, X( KBEG ), 1, CZERO, X( 2*NXFRM+1 ), 1 );
            zgerc(IXFRM, N, -DCMPLX( FACTOR ), X( KBEG ), 1, X( 2*NXFRM+1 ), 1, A( KBEG, 1 ), LDA );

         }

         if ( ITYPE >= 2 && ITYPE <= 4 ) {

            // Apply H(k)* (or H(k)') on the right of A

            if ( ITYPE == 4 ) {
               zlacgv(IXFRM, X( KBEG ), 1 );
            }

            zgemv('N', M, IXFRM, CONE, A( 1, KBEG ), LDA, X( KBEG ), 1, CZERO, X( 2*NXFRM+1 ), 1 );
            zgerc(M, IXFRM, -DCMPLX( FACTOR ), X( 2*NXFRM+1 ), 1, X( KBEG ), 1, A( 1, KBEG ), LDA );

         }
      } // 30

      X( 1 ) = ZLARND( 3, ISEED );
      XABS = ( X( 1 ) ).abs();
      if ( XABS != ZERO ) {
         CSIGN = X( 1 ) / XABS;
      } else {
         CSIGN = CONE;
      }
      X( 2*NXFRM ) = CSIGN;

      // Scale the matrix A by D.

      if ( ITYPE == 1 || ITYPE == 3 || ITYPE == 4 ) {
         for (IROW = 1; IROW <= M; IROW++) { // 40
            zscal(N, DCONJG( X( NXFRM+IROW ) ), A( IROW, 1 ), LDA );
         } // 40
      }

      if ( ITYPE == 2 || ITYPE == 3 ) {
         for (JCOL = 1; JCOL <= N; JCOL++) { // 50
            zscal(M, X( NXFRM+JCOL ), A( 1, JCOL ), 1 );
         } // 50
      }

      if ( ITYPE == 4 ) {
         for (JCOL = 1; JCOL <= N; JCOL++) { // 60
            zscal(M, DCONJG( X( NXFRM+JCOL ) ), A( 1, JCOL ), 1 );
         } // 60
      }
      return;
      }
